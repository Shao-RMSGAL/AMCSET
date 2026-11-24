// Copyright 2025, Texas A&M University
//
// This file is part of AMCSET.
//
// AMCSET is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.

/*! \file winsample.cpp
 *  \brief A minimal sampling profiler for Windows (no administrator rights,
 *  no Visual Studio profiler needed).
 *
 * Starts a program, samples the instruction pointer of each of its threads
 * about every millisecond (SuspendThread / GetThreadContext), and prints the
 * hottest functions (innermost inlined function, from the PDB through
 * DbgHelp), the hottest source lines and the outer (non-inlined) functions.
 * Build from an "x64 Native Tools Command Prompt":
 *
 *     cl /O2 /EHsc /std:c++20 tools\perf\winsample.cpp dbghelp.lib winmm.lib
 *     winsample.exe [--top N] -- build\Release\amcset.exe trim case\TRIM.IN -n 100 -j 1 -q
 *
 * The program needs debug information (Release builds write a PDB).
 * Samples in system DLLs (waits, the C runtime's exp/log) are reported by
 * module.
 */

#define NOMINMAX
#include <windows.h>
#include <dbghelp.h>
#include <psapi.h>
#include <tlhelp32.h>

#include <algorithm>
#include <cstring>
#include <cstdio>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

struct Module {
  std::string path;
  DWORD64 base = 0;
  DWORD size = 0;
};

static std::string module_name(const std::vector<Module> &mods, DWORD64 a) {
  for (const auto &m : mods) {
    if (a >= m.base && a < m.base + m.size) {
      const auto p = m.path.find_last_of("\\/");
      return p == std::string::npos ? m.path : m.path.substr(p + 1);
    }
  }
  return "?";
}

int main(int argc, char **argv) {
  int top = 40;
  int first = 1;
  for (; first < argc; ++first) {
    const std::string a = argv[first];
    if (a == "--top" && first + 1 < argc) top = std::atoi(argv[++first]);
    else if (a == "--") { ++first; break; }
  }
  if (first >= argc) {
    std::fprintf(stderr, "usage: winsample [--top N] -- program args...\n");
    return 2;
  }
  std::string cmd;
  for (int i = first; i < argc; ++i) {
    if (!cmd.empty()) cmd += ' ';
    const std::string a = argv[i];
    cmd += a.find(' ') != std::string::npos ? "\"" + a + "\"" : a;
  }
  STARTUPINFOA si{};
  si.cb = sizeof si;
  PROCESS_INFORMATION pi{};
  if (!CreateProcessA(nullptr, cmd.data(), nullptr, nullptr, FALSE, 0, nullptr, nullptr, &si, &pi)) {
    std::fprintf(stderr, "CreateProcess failed (%lu)\n", GetLastError());
    return 1;
  }
  timeBeginPeriod(1);
  // Stacks of the samples (innermost first, up to kDepth frames).
  constexpr int kDepth = 6;
  std::map<std::vector<DWORD64>, unsigned> hits;
  SymSetOptions(SYMOPT_LOAD_LINES | SYMOPT_UNDNAME | SYMOPT_DEFERRED_LOADS);
  SymInitialize(pi.hProcess, nullptr, FALSE);
  std::map<DWORD64, bool> loaded;
  std::map<DWORD, HANDLE> threads;
  std::vector<Module> mods;
  unsigned long long samples = 0;
  for (int tick = 0;; ++tick) {
    if (WaitForSingleObject(pi.hProcess, 1) == WAIT_OBJECT_0) break;
    if (tick % 50 == 0) {
      HANDLE snap = CreateToolhelp32Snapshot(TH32CS_SNAPTHREAD, 0);
      THREADENTRY32 te{};
      te.dwSize = sizeof te;
      for (BOOL ok = Thread32First(snap, &te); ok; ok = Thread32Next(snap, &te)) {
        if (te.th32OwnerProcessID != pi.dwProcessId || threads.count(te.th32ThreadID)) continue;
        HANDLE h = OpenThread(THREAD_SUSPEND_RESUME | THREAD_GET_CONTEXT | THREAD_QUERY_INFORMATION,
                              FALSE, te.th32ThreadID);
        if (h) threads[te.th32ThreadID] = h;
      }
      CloseHandle(snap);
      HMODULE hm[1024];
      DWORD need = 0;
      if (EnumProcessModules(pi.hProcess, hm, sizeof hm, &need)) {
        mods.clear();
        for (DWORD i = 0; i < need / sizeof(HMODULE); ++i) {
          char path[MAX_PATH];
          MODULEINFO mi{};
          GetModuleFileNameExA(pi.hProcess, hm[i], path, MAX_PATH);
          GetModuleInformation(pi.hProcess, hm[i], &mi, sizeof mi);
          mods.push_back({path, reinterpret_cast<DWORD64>(mi.lpBaseOfDll), mi.SizeOfImage});
          if (!loaded[mods.back().base]) {
            SymLoadModuleEx(pi.hProcess, nullptr, path, nullptr, mods.back().base, mi.SizeOfImage, nullptr, 0);
            loaded[mods.back().base] = true;
          }
        }
      }
    }
    for (auto &[id, h] : threads) {
      if (SuspendThread(h) == static_cast<DWORD>(-1)) continue;
      CONTEXT c{};
      c.ContextFlags = CONTEXT_FULL;
      if (GetThreadContext(h, &c)) {
        std::vector<DWORD64> stack;
        STACKFRAME64 f{};
        f.AddrPC.Offset = c.Rip;
        f.AddrPC.Mode = AddrModeFlat;
        f.AddrFrame.Offset = c.Rbp;
        f.AddrFrame.Mode = AddrModeFlat;
        f.AddrStack.Offset = c.Rsp;
        f.AddrStack.Mode = AddrModeFlat;
        for (int d = 0; d < kDepth; ++d) {
          if (!StackWalk64(IMAGE_FILE_MACHINE_AMD64, pi.hProcess, h, &f, &c, nullptr,
                           SymFunctionTableAccess64, SymGetModuleBase64, nullptr) ||
              f.AddrPC.Offset == 0) {
            break;
          }
          // Return addresses point after the call: step back into it.
          stack.push_back(d == 0 ? f.AddrPC.Offset : f.AddrPC.Offset - 1);
        }
        if (stack.empty()) stack.push_back(c.Rip);
        ++hits[stack];
      }
      ResumeThread(h);
    }
  }
  timeEndPeriod(1);

  HANDLE sym = pi.hProcess;
  char buf[sizeof(SYMBOL_INFO) + 1024];
  // (innermost function, source line, outer function) of an address.
  struct Where {
    std::string inner, line, outer;
  };
  std::unordered_map<DWORD64, Where> cache;
  const auto where = [&](DWORD64 addr) -> const Where & {
    auto it = cache.find(addr);
    if (it != cache.end()) return it->second;
    auto *si2 = reinterpret_cast<SYMBOL_INFO *>(buf);
    std::memset(buf, 0, sizeof buf);
    si2->SizeOfStruct = sizeof(SYMBOL_INFO);
    si2->MaxNameLen = 1000;
    DWORD64 disp = 0;
    std::string outer;
    if (SymFromAddr(sym, addr, &disp, si2)) outer = si2->Name;
    else outer = "[" + module_name(mods, addr) + "]";
    std::string inner = outer;
    std::string line;
    const DWORD inl = SymAddrIncludeInlineTrace(sym, addr);
    DWORD ctx = 0, frame = 0;
    if (inl > 0 && SymQueryInlineTrace(sym, addr, 0, addr, addr, &ctx, &frame)) {
      DWORD64 d2 = 0;
      if (SymFromInlineContext(sym, addr, ctx, &d2, si2)) inner = si2->Name;
      IMAGEHLP_LINE64 li{};
      li.SizeOfStruct = sizeof li;
      DWORD d3 = 0;
      if (SymGetLineFromInlineContext(sym, addr, ctx, 0, &d3, &li)) {
        std::string f = li.FileName;
        const auto p = f.find_last_of("\\/");
        line = (p == std::string::npos ? f : f.substr(p + 1)) + ":" + std::to_string(li.LineNumber);
      }
    } else {
      IMAGEHLP_LINE64 li{};
      li.SizeOfStruct = sizeof li;
      DWORD d3 = 0;
      if (SymGetLineFromAddr64(sym, addr, &d3, &li)) {
        std::string f = li.FileName;
        const auto p = f.find_last_of("\\/");
        line = (p == std::string::npos ? f : f.substr(p + 1)) + ":" + std::to_string(li.LineNumber);
      }
    }
    if (line.empty()) line = "[" + module_name(mods, addr) + "]";
    return cache[addr] = Where{inner, line, outer};
  };
  const auto idle = [](const std::string &f) {
    return f.rfind("NtWait", 0) == 0 || f.rfind("ZwWait", 0) == 0 || f.rfind("ZwDelay", 0) == 0 ||
           f.rfind("NtDelay", 0) == 0 || f.rfind("NtRemoveIo", 0) == 0 || f.rfind("ZwRemoveIo", 0) == 0 ||
           f.rfind("NtWaitForWorkViaWorkerFactory", 0) == 0;
  };
  std::map<std::string, unsigned long long> by_func, by_line, by_outer, by_caller, inclusive;
  for (const auto &[stack, n] : hits) {
    const Where &w0 = where(stack[0]);
    if (idle(w0.inner)) continue;
    samples += n;
    by_func[w0.inner] += n;
    by_outer[w0.outer] += n;
    by_line[w0.line + "  " + w0.inner] += n;
    // Leaf functions outside the program (C runtime math, allocation):
    // charge them to the first caller in amcset code.
    if (w0.line.front() == '[') {
      for (std::size_t d = 1; d < stack.size(); ++d) {
        const Where &w = where(stack[d]);
        if (w.line.front() != '[') {
          by_caller[w0.inner + "  <-  " + w.line + " " + w.inner] += n;
          break;
        }
      }
    }
    std::map<std::string, bool> seen;
    for (const auto a : stack) {
      const Where &w = where(a);
      if (!seen[w.outer]) inclusive[w.outer] += n;
      seen[w.outer] = true;
    }
  }
  const auto print = [&](const char *title, const std::map<std::string, unsigned long long> &m) {
    std::vector<std::pair<unsigned long long, std::string>> v;
    for (const auto &[k, n] : m) v.push_back({n, k});
    std::sort(v.rbegin(), v.rend());
    std::printf("\n== %s (%llu samples) ==\n", title, samples);
    for (int i = 0; i < top && i < static_cast<int>(v.size()); ++i) {
      std::printf("%6.2f%%  %s\n", 100.0 * static_cast<double>(v[i].first) / static_cast<double>(samples),
                  v[i].second.substr(0, 150).c_str());
    }
  };
  print("innermost function (busy samples)", by_func);
  print("source line", by_line);
  print("outer function", by_outer);
  print("library leaf <- caller", by_caller);
  print("inclusive outer function (stack depth 6)", inclusive);
  DWORD code = 0;
  GetExitCodeProcess(pi.hProcess, &code);
  return static_cast<int>(code);
}

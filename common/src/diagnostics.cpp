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

/*! \file diagnostics.cpp
 *  \brief Diagnostic log, crash reports and breadcrumbs (see diagnostics.h).
 *
 * The crash path avoids the heap and locks where it can: the report is
 * built in a static buffer with hand-written number formatting, the log
 * ring and the breadcrumbs are fixed arrays read without locks (a report
 * may show a line being written at that moment; it never blocks). On
 * Windows the report is written by a fresh thread, so that a stack overflow
 * still gets one.
 */

#include "diagnostics.h"

#include "config.h"

#include <absl/log/log_entry.h>
#include <absl/log/log_sink.h>
#include <absl/log/log_sink_registry.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
// dbghelp.h needs windows.h first.
#include <dbghelp.h>
#include <io.h>
#pragma comment(lib, "dbghelp.lib")
#else
#include <absl/debugging/stacktrace.h>
#include <absl/debugging/symbolize.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#endif

namespace amcset {
namespace diag {

namespace {

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

std::atomic<int> g_level{static_cast<int>(Level::warning)};
std::mutex g_mutex; // log file, options
std::FILE *g_file = nullptr;
bool g_echo = false;
std::filesystem::path g_log_path;

constexpr std::size_t kPathMax = 1024;
char g_program[64] = "amcset";
char g_crash_dir[kPathMax] = "";
std::filesystem::path g_crash_dir_path;

// The last log lines (crash reports).
constexpr std::size_t kRingLines = 48;
constexpr std::size_t kRingWidth = 240;
char g_ring[kRingLines][kRingWidth];
std::atomic<std::uint64_t> g_ring_next{0};

// Run context.
constexpr std::size_t kContextMax = 2048;
char g_context[kContextMax] = "";
std::mutex g_context_mutex;

// Breadcrumbs.
constexpr int kSlots = 256;
constexpr std::uint64_t kNone = ~std::uint64_t{0};
struct Slot {
  std::atomic<bool> used{false};
  char what[48] = {};
  std::atomic<std::uint64_t> thread{0};
  std::atomic<std::uint64_t> particle{kNone};
  std::atomic<std::uint64_t> block_first{0};
  std::atomic<std::uint32_t> seed{0};
};
Slot g_slots[kSlots];

std::atomic<bool> g_crashing{false};
std::atomic<bool> g_handlers_installed{false};
// Windows: write a minidump with each crash report (environment variable
// AMCSET_MINIDUMP=0 turns it off, e.g. for the tests).
bool g_minidump = true;

std::uint64_t current_thread_number() noexcept {
#if defined(_WIN32)
  return static_cast<std::uint64_t>(GetCurrentThreadId());
#else
  return static_cast<std::uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id()));
#endif
}

// ---------------------------------------------------------------------------
// Allocation-free text building for the crash path
// ---------------------------------------------------------------------------

struct Text {
  char *buf;
  std::size_t cap;
  std::size_t len = 0;
  void add(const char *s) noexcept {
    while (s && *s && len + 1 < cap) buf[len++] = *s++;
    buf[len] = 0;
  }
  void add(const char *s, std::size_t n) noexcept {
    for (std::size_t i = 0; i < n && s[i] && len + 1 < cap; ++i) buf[len++] = s[i];
    buf[len] = 0;
  }
  void num(std::uint64_t v) noexcept {
    char t[24];
    int n = 0;
    do {
      t[n++] = static_cast<char>('0' + v % 10);
      v /= 10;
    } while (v != 0 && n < 24);
    while (n > 0 && len + 1 < cap) buf[len++] = t[--n];
    buf[len] = 0;
  }
  void hex(std::uint64_t v) noexcept {
    add("0x");
    char t[17];
    int n = 0;
    do {
      t[n++] = "0123456789abcdef"[v & 15];
      v >>= 4;
    } while (v != 0 && n < 16);
    while (n > 0 && len + 1 < cap) buf[len++] = t[--n];
    buf[len] = 0;
  }
  void two(unsigned v) noexcept {
    if (v < 10) add("0");
    num(v);
  }
};

//! UTC date and time of \p t (seconds since 1970), without gmtime (not
//! async-signal-safe). Howard Hinnant's civil-from-days algorithm.
void utc_stamp(Text &out, std::int64_t t, bool compact) noexcept {
  std::int64_t days = t / 86400;
  std::int64_t sec = t % 86400;
  if (sec < 0) {
    sec += 86400;
    --days;
  }
  days += 719468;
  const std::int64_t era = (days >= 0 ? days : days - 146096) / 146097;
  const auto doe = static_cast<unsigned>(days - era * 146097);
  const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const std::int64_t y = static_cast<std::int64_t>(yoe) + era * 400;
  const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const unsigned mp = (5 * doy + 2) / 153;
  const unsigned d = doy - (153 * mp + 2) / 5 + 1;
  const unsigned m = mp < 10 ? mp + 3 : mp - 9;
  out.num(static_cast<std::uint64_t>(y + (m <= 2 ? 1 : 0)));
  if (!compact) out.add("-");
  out.two(m);
  if (!compact) out.add("-");
  out.two(d);
  out.add(compact ? "-" : " ");
  out.two(static_cast<unsigned>(sec / 3600));
  if (!compact) out.add(":");
  out.two(static_cast<unsigned>(sec / 60 % 60));
  if (!compact) out.add(":");
  out.two(static_cast<unsigned>(sec % 60));
  if (!compact) out.add(" UTC");
}

std::uint64_t process_id() noexcept {
#if defined(_WIN32)
  return GetCurrentProcessId();
#else
  return static_cast<std::uint64_t>(getpid());
#endif
}

char g_report[64 * 1024];

//! The parts of a report that do not depend on the platform.
void report_header(Text &r, const char *what) noexcept {
  r.add("==== AMCSET crash report ====\nProgram:  ");
  r.add(g_program);
  r.add(" (AMCSET ");
  r.add(AMCSET_VERSION);
  r.add(")\nTime:     ");
  utc_stamp(r, static_cast<std::int64_t>(std::time(nullptr)), false);
  r.add("\nProcess:  ");
  r.num(process_id());
  r.add("  crashing thread ");
  r.num(current_thread_number());
  r.add("\nProblem:  ");
  r.add(what);
  r.add("\n");
}

void report_context(Text &r) noexcept {
  r.add("\n-- Calculation --\n");
  if (g_context[0] != 0) {
    r.add(g_context, kContextMax);
    r.add("\n");
  } else {
    r.add("(none recorded)\n");
  }
  r.add("\n-- Threads at work --\n");
  bool any = false;
  for (auto &s : g_slots) {
    if (!s.used.load(std::memory_order_relaxed)) continue;
    any = true;
    r.add("thread ");
    r.num(s.thread.load(std::memory_order_relaxed));
    r.add(" (");
    r.add(s.what, sizeof s.what);
    r.add("): ");
    const std::uint64_t p = s.particle.load(std::memory_order_relaxed);
    if (p == kNone) {
      r.add("idle\n");
      continue;
    }
    const std::uint64_t first = s.block_first.load(std::memory_order_relaxed);
    r.add("particle ");
    r.num(p);
    r.add(" (0-based), seed ");
    r.num(s.seed.load(std::memory_order_relaxed));
    r.add("\n    rerun exactly with: --set first_ion=");
    r.num(first);
    r.add(" -n ");
    r.num(p - first + 1);
    r.add(" -j 1 (same scenario and seed)\n");
  }
  if (!any) r.add("(no worker registered)\n");
  r.add("\n-- Last log lines --\n");
  const std::uint64_t next = g_ring_next.load(std::memory_order_relaxed);
  const std::uint64_t first = next > kRingLines ? next - kRingLines : 0;
  if (next == 0) r.add("(log empty; run with --log FILE or the GUI's diagnostic log for more)\n");
  for (std::uint64_t i = first; i < next; ++i) {
    r.add(g_ring[i % kRingLines], kRingWidth);
    r.add("\n");
  }
}

#if defined(_WIN32)

// ---------------------------------------------------------------------------
// Windows
// ---------------------------------------------------------------------------

struct CrashInfo {
  const char *what = nullptr;
  EXCEPTION_POINTERS *pointers = nullptr;
  DWORD thread_id = 0;
  CONTEXT context{};
  bool have_context = false;
};

const char *exception_name(DWORD code) noexcept {
  switch (code) {
  case EXCEPTION_ACCESS_VIOLATION: return "access violation (invalid memory access)";
  case EXCEPTION_STACK_OVERFLOW: return "stack overflow";
  case EXCEPTION_INT_DIVIDE_BY_ZERO: return "integer division by zero";
  case EXCEPTION_ILLEGAL_INSTRUCTION: return "illegal instruction";
  case EXCEPTION_ARRAY_BOUNDS_EXCEEDED: return "array bounds exceeded";
  case EXCEPTION_IN_PAGE_ERROR: return "in-page error (file or memory mapping failed)";
  case EXCEPTION_PRIV_INSTRUCTION: return "privileged instruction";
  case EXCEPTION_FLT_DIVIDE_BY_ZERO: return "floating-point division by zero";
  case EXCEPTION_FLT_INVALID_OPERATION: return "invalid floating-point operation";
  case 0xE06D7363: return "unhandled C++ exception";
  default: return "unhandled exception";
  }
}

void open_report_file(Text &name, const char *ext) noexcept {
  name.add(g_crash_dir);
  name.add("\\");
  name.add(g_program);
  name.add("-crash-");
  utc_stamp(name, static_cast<std::int64_t>(std::time(nullptr)), true);
  name.add("-");
  name.num(process_id());
  name.add(ext);
}

void append_stack(Text &r, HANDLE thread, CONTEXT ctx) noexcept {
  HANDLE process = GetCurrentProcess();
  SymSetOptions(SYMOPT_LOAD_LINES | SYMOPT_UNDNAME | SYMOPT_DEFERRED_LOADS);
  SymInitialize(process, nullptr, TRUE);
  STACKFRAME64 f{};
  f.AddrPC.Offset = ctx.Rip;
  f.AddrPC.Mode = AddrModeFlat;
  f.AddrFrame.Offset = ctx.Rbp;
  f.AddrFrame.Mode = AddrModeFlat;
  f.AddrStack.Offset = ctx.Rsp;
  f.AddrStack.Mode = AddrModeFlat;
  alignas(SYMBOL_INFO) char symbuf[sizeof(SYMBOL_INFO) + 512];
  for (int depth = 0; depth < 64; ++depth) {
    if (!StackWalk64(IMAGE_FILE_MACHINE_AMD64, process, thread, &f, &ctx, nullptr,
                     SymFunctionTableAccess64, SymGetModuleBase64, nullptr) ||
        f.AddrPC.Offset == 0) {
      break;
    }
    const DWORD64 pc = f.AddrPC.Offset;
    r.add("  #");
    r.num(static_cast<std::uint64_t>(depth));
    r.add(" ");
    r.hex(pc);
    auto *sym = reinterpret_cast<SYMBOL_INFO *>(symbuf);
    std::memset(symbuf, 0, sizeof symbuf);
    sym->SizeOfStruct = sizeof(SYMBOL_INFO);
    sym->MaxNameLen = 500;
    DWORD64 disp = 0;
    if (SymFromAddr(process, pc, &disp, sym)) {
      r.add(" ");
      r.add(sym->Name);
    }
    IMAGEHLP_LINE64 line{};
    line.SizeOfStruct = sizeof line;
    DWORD ldisp = 0;
    if (SymGetLineFromAddr64(process, depth == 0 ? pc : pc - 1, &ldisp, &line)) {
      r.add("  ");
      r.add(line.FileName);
      r.add(":");
      r.num(line.LineNumber);
    } else {
      IMAGEHLP_MODULE64 mod{};
      mod.SizeOfStruct = sizeof mod;
      if (SymGetModuleInfo64(process, pc, &mod)) {
        r.add("  [");
        r.add(mod.ModuleName);
        r.add("+");
        r.hex(pc - mod.BaseOfImage);
        r.add("]");
      }
    }
    r.add("\n");
  }
}

DWORD WINAPI report_thread(LPVOID param) {
  auto *info = static_cast<CrashInfo *>(param);
  Text r{g_report, sizeof g_report};
  report_header(r, info->what);
  if (info->pointers != nullptr && info->pointers->ExceptionRecord != nullptr) {
    const auto *er = info->pointers->ExceptionRecord;
    r.add("Code:     ");
    r.hex(er->ExceptionCode);
    r.add(" at ");
    r.hex(reinterpret_cast<std::uint64_t>(er->ExceptionAddress));
    if (er->ExceptionCode == EXCEPTION_ACCESS_VIOLATION && er->NumberParameters >= 2) {
      r.add(er->ExceptionInformation[0] == 0 ? " (reading " : " (writing ");
      r.hex(er->ExceptionInformation[1]);
      r.add(")");
    }
    r.add("\n");
  }
  r.add("\n-- Stack of the crashing thread --\n");
  if (info->have_context) {
    HANDLE thread = OpenThread(THREAD_ALL_ACCESS, FALSE, info->thread_id);
    append_stack(r, thread != nullptr ? thread : GetCurrentThread(), info->context);
    if (thread != nullptr) CloseHandle(thread);
  } else {
    r.add("(not available)\n");
  }
  report_context(r);
  // Files: report and minidump.
  char namebuf[kPathMax + 128];
  Text name{namebuf, sizeof namebuf};
  open_report_file(name, ".txt");
  HANDLE f = CreateFileA(namebuf, GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
  r.add("\nReport file: ");
  r.add(namebuf);
  r.add("\n");
  char dumpbuf[kPathMax + 128];
  Text dump{dumpbuf, sizeof dumpbuf};
  open_report_file(dump, ".dmp");
  r.add("Minidump:    ");
  r.add(g_minidump ? dumpbuf : "(not written: AMCSET_MINIDUMP=0)");
  r.add("\n");
  if (f != INVALID_HANDLE_VALUE) {
    DWORD written = 0;
    WriteFile(f, r.buf, static_cast<DWORD>(r.len), &written, nullptr);
    CloseHandle(f);
  }
  HANDLE d = g_minidump ? CreateFileA(dumpbuf, GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS,
                                      FILE_ATTRIBUTE_NORMAL, nullptr)
                        : INVALID_HANDLE_VALUE;
  if (d != INVALID_HANDLE_VALUE) {
    MINIDUMP_EXCEPTION_INFORMATION mei{};
    mei.ThreadId = info->thread_id;
    mei.ExceptionPointers = info->pointers;
    mei.ClientPointers = FALSE;
    const auto type = static_cast<MINIDUMP_TYPE>(MiniDumpWithIndirectlyReferencedMemory |
                                                 MiniDumpWithDataSegs | MiniDumpWithThreadInfo);
    MiniDumpWriteDump(GetCurrentProcess(), GetCurrentProcessId(), d, type,
                      info->pointers != nullptr ? &mei : nullptr, nullptr, nullptr);
    CloseHandle(d);
  }
  // stderr last: the report is on disk even if the console is gone.
  DWORD n = 0;
  HANDLE err = GetStdHandle(STD_ERROR_HANDLE);
  if (err != nullptr && err != INVALID_HANDLE_VALUE) WriteFile(err, r.buf, static_cast<DWORD>(r.len), &n, nullptr);
  {
    // The log file as well (best effort, without its mutex).
    if (g_file != nullptr) {
      std::fwrite(r.buf, 1, r.len, g_file);
      std::fflush(g_file);
    }
  }
  return 0;
}

//! Write a report from a fresh thread (works after a stack overflow too).
void crash_report(const char *what, EXCEPTION_POINTERS *pointers, const CONTEXT *ctx) noexcept {
  if (g_crashing.exchange(true)) return; // one report per process
  static CrashInfo info;
  info.what = what;
  info.pointers = pointers;
  info.thread_id = GetCurrentThreadId();
  if (ctx != nullptr) {
    info.context = *ctx;
    info.have_context = true;
  }
  HANDLE t = CreateThread(nullptr, 1 << 20, report_thread, &info, 0, nullptr);
  if (t != nullptr) {
    WaitForSingleObject(t, 60000);
    CloseHandle(t);
  } else {
    report_thread(&info);
  }
}

LONG WINAPI unhandled_filter(EXCEPTION_POINTERS *p) {
  crash_report(exception_name(p->ExceptionRecord->ExceptionCode), p, p->ContextRecord);
  // Terminate without Windows Error Reporting's dialog (it would hang
  // unattended runs), with exit code 3 as after abort(); the report names
  // the exception code.
  TerminateProcess(GetCurrentProcess(), 3);
  return EXCEPTION_EXECUTE_HANDLER;
}

void report_here(const char *what) noexcept {
  CONTEXT ctx{};
  RtlCaptureContext(&ctx);
  crash_report(what, nullptr, &ctx);
}

#else

// ---------------------------------------------------------------------------
// POSIX
// ---------------------------------------------------------------------------

void write_all(int fd, const char *p, std::size_t n) noexcept {
  while (n > 0) {
    const ssize_t w = ::write(fd, p, n);
    if (w <= 0) return;
    p += w;
    n -= static_cast<std::size_t>(w);
  }
}

void append_stack(Text &r, int skip) noexcept {
  void *frames[64];
  const int n = absl::GetStackTrace(frames, 64, skip);
  char sym[512];
  for (int i = 0; i < n; ++i) {
    r.add("  #");
    r.num(static_cast<std::uint64_t>(i));
    r.add(" ");
    r.hex(reinterpret_cast<std::uint64_t>(frames[i]));
    // Return addresses point after the call.
    void *pc = i == 0 ? frames[i] : static_cast<char *>(frames[i]) - 1;
    if (absl::Symbolize(pc, sym, sizeof sym)) {
      r.add(" ");
      r.add(sym);
    }
    r.add("\n");
  }
}

void crash_report(const char *what, int skip) noexcept {
  if (g_crashing.exchange(true)) return;
  Text r{g_report, sizeof g_report};
  report_header(r, what);
  r.add("\n-- Stack of the crashing thread --\n");
  append_stack(r, skip);
  report_context(r);
  char namebuf[kPathMax + 128];
  Text name{namebuf, sizeof namebuf};
  name.add(g_crash_dir);
  name.add("/");
  name.add(g_program);
  name.add("-crash-");
  utc_stamp(name, static_cast<std::int64_t>(std::time(nullptr)), true);
  name.add("-");
  name.num(process_id());
  name.add(".txt");
  r.add("\nReport file: ");
  r.add(namebuf);
  r.add("\n");
  const int fd = ::open(namebuf, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fd >= 0) {
    write_all(fd, r.buf, r.len);
    ::close(fd);
  }
  write_all(2, r.buf, r.len);
  if (g_file != nullptr) write_all(fileno(g_file), r.buf, r.len);
}

const char *signal_name(int sig) noexcept {
  switch (sig) {
  case SIGSEGV: return "SIGSEGV (invalid memory access)";
  case SIGBUS: return "SIGBUS (bus error)";
  case SIGFPE: return "SIGFPE (arithmetic error)";
  case SIGILL: return "SIGILL (illegal instruction)";
  case SIGABRT: return "SIGABRT (abort)";
  default: return "fatal signal";
  }
}

void signal_handler(int sig, siginfo_t *, void *) {
  crash_report(signal_name(sig), 1);
  // SA_RESETHAND restored the default action: terminate (and dump core).
  ::raise(sig);
}

alignas(16) char g_altstack[64 * 1024];

void report_here(const char *what) noexcept { crash_report(what, 1); }

#endif

void terminate_handler() {
  // Name the exception that got here, if any.
  static char what[512];
  Text t{what, sizeof what};
  t.add("std::terminate");
  if (auto e = std::current_exception()) {
    try {
      std::rethrow_exception(e);
    } catch (const std::exception &x) {
      t.add(": uncaught exception: ");
      t.add(x.what());
    } catch (...) {
      t.add(": uncaught exception of unknown type");
    }
  } else {
    t.add(" (no active exception)");
  }
  report_here(what);
  std::_Exit(3);
}

#if defined(_WIN32)
void purecall_handler() {
  report_here("pure virtual function call");
  std::_Exit(3);
}

void invalid_parameter_handler(const wchar_t *, const wchar_t *, const wchar_t *, unsigned, uintptr_t) {
  report_here("invalid parameter passed to a C runtime function");
  std::_Exit(3);
}

void abort_handler(int) {
  report_here("abort() called");
  std::_Exit(3);
}
#endif

// ---------------------------------------------------------------------------
// abseil messages into the log
// ---------------------------------------------------------------------------

class AbslSink final : public absl::LogSink {
public:
  void Send(const absl::LogEntry &e) override {
    Level l = Level::info;
    switch (e.log_severity()) {
    case absl::LogSeverity::kInfo: l = Level::info; break;
    case absl::LogSeverity::kWarning: l = Level::warning; break;
    default: l = Level::error; break;
    }
    if (!enabled(l)) return;
    std::string msg = "[abseil ";
    msg += std::string(e.source_basename());
    msg += ":" + std::to_string(e.source_line()) + "] ";
    msg += std::string(e.text_message());
    write(l, msg);
  }
};

void set_crash_dir(const std::filesystem::path &dir) {
  std::error_code ec;
  std::filesystem::create_directories(dir, ec);
  g_crash_dir_path = dir;
  const std::string s = dir.string();
  std::snprintf(g_crash_dir, sizeof g_crash_dir, "%s", s.c_str());
}

} // namespace

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

const char *to_string(Level level) noexcept {
  switch (level) {
  case Level::debug: return "debug";
  case Level::info: return "info";
  case Level::warning: return "warning";
  case Level::error: return "error";
  case Level::off: return "off";
  }
  return "off";
}

Level parse_level(std::string_view text) {
  if (text == "debug") return Level::debug;
  if (text == "info") return Level::info;
  if (text == "warning" || text == "warn") return Level::warning;
  if (text == "error") return Level::error;
  if (text == "off" || text == "none") return Level::off;
  throw std::invalid_argument("unknown log level '" + std::string(text) +
                              "' (expected debug, info, warning, error or off)");
}

void set_level(Level level) noexcept { g_level.store(static_cast<int>(level), std::memory_order_relaxed); }

Level level() noexcept { return static_cast<Level>(g_level.load(std::memory_order_relaxed)); }

bool enabled(Level l) noexcept {
  return l != Level::off && static_cast<int>(l) >= g_level.load(std::memory_order_relaxed);
}

void write(Level l, std::string_view message) {
  if (!enabled(l)) return;
  // "2026-11-23 18:04:05.123 UTC [info] [thread 1234] message"
  char stamp[64];
  Text ts{stamp, sizeof stamp};
  const auto now = std::chrono::system_clock::now();
  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
  utc_stamp(ts, static_cast<std::int64_t>(ms / 1000), false);
  char frac[8];
  std::snprintf(frac, sizeof frac, ".%03d", static_cast<int>(ms % 1000));
  std::string line = std::string(stamp, ts.len >= 4 ? ts.len - 4 : ts.len); // drop " UTC"
  line += frac;
  line += " [";
  line += to_string(l);
  line += "] [thread " + std::to_string(current_thread_number()) + "] ";
  line.append(message.data(), message.size());
  while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
  {
    const std::uint64_t k = g_ring_next.fetch_add(1, std::memory_order_acq_rel);
    std::snprintf(g_ring[k % kRingLines], kRingWidth, "%s", line.c_str());
  }
  line += '\n';
  std::lock_guard<std::mutex> lock(g_mutex);
  if (g_file != nullptr) {
    std::fwrite(line.data(), 1, line.size(), g_file);
    // Flushed every line: the log must survive a crash. Lines are rare
    // (run-level events), so this costs nothing measurable.
    std::fflush(g_file);
  }
  if (g_echo) {
    std::fwrite(line.data(), 1, line.size(), stderr);
  }
}

std::filesystem::path log_file() {
  std::lock_guard<std::mutex> lock(g_mutex);
  return g_log_path;
}

std::filesystem::path crash_directory() {
  std::lock_guard<std::mutex> lock(g_mutex);
  return g_crash_dir_path;
}

void set_context(std::string_view text) {
  std::lock_guard<std::mutex> lock(g_context_mutex);
  const std::size_t n = std::min(text.size(), kContextMax - 1);
  std::memcpy(g_context, text.data(), n);
  g_context[n] = 0;
}

void install_crash_handlers(const char *argv0) {
  if (g_handlers_installed.exchange(true)) return;
  {
    std::lock_guard<std::mutex> lock(g_mutex);
    if (g_crash_dir[0] == 0) {
      std::error_code ec;
      auto dir = std::filesystem::temp_directory_path(ec);
      if (ec) dir = ".";
      set_crash_dir(dir / "amcset-crashes");
    }
  }
  std::set_terminate(terminate_handler);
  if (const char *e = std::getenv("AMCSET_MINIDUMP"); e != nullptr && std::strcmp(e, "0") == 0) {
    g_minidump = false;
  }
#if defined(_WIN32)
  (void)argv0;
  SetUnhandledExceptionFilter(unhandled_filter);
  _set_purecall_handler(purecall_handler);
  _set_invalid_parameter_handler(invalid_parameter_handler);
  // abort(): our report instead of the CRT's message box / Windows Error
  // Reporting (which would hang unattended runs).
  _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
  std::signal(SIGABRT, abort_handler);
#else
  absl::InitializeSymbolizer(argv0 != nullptr ? argv0 : "amcset");
  stack_t ss{};
  ss.ss_sp = g_altstack;
  ss.ss_size = sizeof g_altstack;
  sigaltstack(&ss, nullptr);
  struct sigaction sa{};
  sa.sa_sigaction = signal_handler;
  sa.sa_flags = SA_SIGINFO | SA_ONSTACK | SA_RESETHAND;
  sigemptyset(&sa.sa_mask);
  for (int sig : {SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGABRT}) sigaction(sig, &sa, nullptr);
#endif
}

void initialize(const char *argv0, const Options &o) {
  {
    std::lock_guard<std::mutex> lock(g_mutex);
    std::snprintf(g_program, sizeof g_program, "%s", o.program.c_str());
    g_echo = o.echo_to_stderr;
    if (g_file != nullptr) {
      std::fclose(g_file);
      g_file = nullptr;
    }
    g_log_path.clear();
    if (!o.log_file.empty()) {
      std::error_code ec;
      if (o.log_file.has_parent_path()) std::filesystem::create_directories(o.log_file.parent_path(), ec);
#if defined(_WIN32)
      g_file = _wfopen(o.log_file.c_str(), L"ab");
#else
      g_file = std::fopen(o.log_file.c_str(), "ab");
#endif
      if (g_file == nullptr) {
        std::fprintf(stderr, "amcset: cannot open the log file %s; logging to memory only\n",
                     o.log_file.string().c_str());
      } else {
        g_log_path = o.log_file;
      }
    }
    if (!o.crash_directory.empty()) {
      set_crash_dir(o.crash_directory);
    } else if (!g_log_path.empty()) {
      std::error_code ec;
      const auto abs = std::filesystem::absolute(g_log_path, ec);
      set_crash_dir(ec ? g_log_path.parent_path() : abs.parent_path());
    }
  }
  set_level(o.level);
  static AbslSink sink;
  static std::once_flag sink_once;
  std::call_once(sink_once, [] { absl::AddLogSink(&sink); });
  install_crash_handlers(argv0);
}

Activity::Activity(const char *what) noexcept {
  for (int i = 0; i < kSlots; ++i) {
    bool expected = false;
    if (g_slots[i].used.compare_exchange_strong(expected, true)) {
      slot_ = i;
      auto &s = g_slots[i];
      std::snprintf(s.what, sizeof s.what, "%s", what != nullptr ? what : "");
      s.thread.store(current_thread_number(), std::memory_order_relaxed);
      s.particle.store(kNone, std::memory_order_relaxed);
      return;
    }
  }
}

Activity::~Activity() {
  if (slot_ >= 0) {
    g_slots[slot_].particle.store(kNone, std::memory_order_relaxed);
    g_slots[slot_].used.store(false, std::memory_order_release);
  }
}

void Activity::set_particle(std::uint64_t index, std::uint64_t block_first, std::uint32_t seed) noexcept {
  if (slot_ < 0) return;
  auto &s = g_slots[slot_];
  s.block_first.store(block_first, std::memory_order_relaxed);
  s.seed.store(seed, std::memory_order_relaxed);
  s.particle.store(index, std::memory_order_relaxed);
}

void Activity::clear_particle() noexcept {
  if (slot_ >= 0) g_slots[slot_].particle.store(kNone, std::memory_order_relaxed);
}

} // namespace diag
} // namespace amcset

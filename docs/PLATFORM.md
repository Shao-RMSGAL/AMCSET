# Platform requirements and how AMCSET meets them

What Windows and Linux expect of a desktop and command-line program, what
AMCSET 1.0 does about each, and how that was checked. Packaging and
signing: `docs/RELEASING.md`; installing: `docs/INSTALL.md`.

## Windows

The list follows Microsoft's "Desktop certification requirements" for
Windows 10 desktop apps (document version 10, 2015; archived on Microsoft
Learn:
https://learn.microsoft.com/en-us/previous-versions/windows/win32/win_cert/certification-requirements-for-windows-desktop-apps),
the requirements the Windows App Certification Kit tests, plus Smart App
Control's signing rule (RELEASING.md §3).

| Requirement | AMCSET 1.0 | Checked |
|---|---|---|
| **1. Compatibility and resiliency**: a compatibility manifest; DPI awareness declared in the manifest; no AppInit_DLLs, no VB6 runtime | `packaging/windows/amcset.manifest` in both executables: `supportedOS` Windows 10/11, `dpiAwareness` Per-Monitor v2 (Qt's own default, which Qt then keeps), `longPathAware`, `activeCodePage` UTF-8 | manifest extracted from the built `amcset.exe` with `mt.exe`; a run with an output path containing "ψ" (not in the old ANSI code page) |
| **3. Windows security features**: /DYNAMICBASE (ASLR), /NXCOMPAT (DEP), no shared writable PE sections (/SafeSEH is for 32-bit code) | MSVC defaults, plus high-entropy ASLR and /GS stack checks; Control Flow Guard is not enabled (not measured for speed yet) | `dumpbin /headers`: "Dynamic base", "NX compatible", "High Entropy Virtual Addresses"; /GS in every object |
| **4. Restart Manager messages**: answer WM_QUERYENDSESSION at once with TRUE, handle critical shutdowns, close within 30 s | GUI: at the session's end it stops a running calculation, keeps a changed scenario in `%LOCALAPPDATA%\AMCSET\AMCSET\recovery\` (reported at the next start) and closes without asking or vetoing (`MainWindow::session_ending`). CLI: console close, log-off and shutdown events stop the run and write the results within the 5 s Windows allows (`cli/src/interrupts.cpp`) | code review; the GUI path was not triggered by a real log-off |
| **5. Clean, reversible installation**: silent install and uninstall, an Add/Remove Programs entry, no forced restart, roll back on failure | MSI (Windows Installer: transactional, `msiexec /qn`, Settings > Apps entry, upgrade code fixed so that new versions replace old ones, no restart); ZIP needs no installation | ZIP and MSIX built and inspected here; the MSI is built by the release workflow (WiX is not installed on the workstation): test one before publishing (RELEASING.md §6) |
| **6. Digitally signed files**: every .exe and .dll, and the installer | `tools/windows/sign.ps1` and the release workflow sign the programs, the Qt DLLs, the MSI and the MSIX (Authenticode, SHA-256, RFC 3161 time stamp) once a certificate is set up | **open: 1.0 is unsigned**; the maintainer chose the SignPath Foundation, to be set up later (RELEASING.md §3.2) |
| **7. No version checks that block installation or launch** | none (the MSIX declares Windows 10 1809 as its minimum, the system's own check) | code review |
| **8. No services or drivers in Safe Mode** | none installed | — |
| **9. User Account Control**: the main process runs asInvoker | `asInvoker` (the linker's manifest part); the per-machine MSI elevates through Windows Installer only | extracted manifest |
| **10. Correct folders**: Program Files by default; per-user data in the user's profile; no auto-start | MSI installs to `C:\Program Files\AMCSET`; settings in the registry (QSettings), logs and recovery copies in `%LOCALAPPDATA%\AMCSET`, CLI crash reports in `%TEMP%`; results only where the user says; nothing starts with Windows | install layout listed from `cmake --install` |
| **11. Multi-user sessions** | nothing machine-wide is written at run time; each user has their own settings | code review |
| **12. Native x64** | x64 only | — |
| **Version information and icon** (shown by Explorer, Task Manager, installers) | `packaging/windows/version.rc.in`: company, product, version 1.0.0, description, copyright; the logo | `(Get-Item amcset.exe).VersionInfo` |
| **Packaging**: MSIX manifest valid, full trust declared | `packaging/windows/AppxManifest.xml.in`: runFullTrust, the GUI in Start, `amcset` and `amcset-gui` as execution aliases (the CLI as a console program) | `makeappx pack` validates the manifest against the schema |
| **Smart App Control / SmartScreen** | signed releases (above); until then docs/INSTALL.md "Running the unsigned Windows release" says what users can do | open, with item 6 |

Not run: the Windows App Certification Kit (`appcert.exe`, installed with
the Windows SDK). It installs the program under test, which needs the MSI
(release workflow) and administrator rights. Run it on a signed MSI before
the first public release (RELEASING.md §6).

## Linux

| Standard | AMCSET 1.0 | Checked |
|---|---|---|
| **Signals** (POSIX): SIGTERM, SIGHUP, SIGINT end a program cleanly; status 128 + n for a program ended by signal n | CLI: the run stops after the ions in progress, partial results are written, then the signal is re-raised with its default action, so the parent sees the signal (a shell reports 143 for SIGTERM, 130 for SIGINT, and a calling script stops); a second signal ends it at once; handlers only set lock-free atomics. GUI: a self-pipe turns them into a quiet quit on the GUI thread (status 128 + n). SIGPIPE and SIGQUIT keep their default actions | `Cli.InterruptsSetTheFlagAndGiveTheSignalsStatus`; in WSL `kill -TERM` during a run: results written, status 143 |
| **Exit status and streams** | 0 success, 1 error, 2 usage, 3 crash, 128+n signal; results and `--help`/`--version` on stdout, messages on stderr | CLI tests |
| **Locale** | numbers in files always use a decimal point: the GUI resets `LC_NUMERIC` to "C" after Qt sets the locale from the environment (as Qt advises); the interface still shows the user's locale | code review |
| **Filesystem Hierarchy Standard** | `/usr/bin/amcset`, `/usr/bin/amcset-gui`, `/usr/share/amcset/`, `/usr/share/doc/amcset/` (with Debian's `copyright` file), `/usr/share/man/man1/` | `.deb` contents listed (WSL) |
| **XDG base directories** | settings `~/.config/AMCSET/`, logs and recovery copies `~/.local/share/AMCSET/AMCSET/` (Qt's standard paths) | — |
| **Desktop Entry specification** | `io.github.Shao_RMSGAL.AMCSET.desktop`; the GUI sets the same desktop file name (Wayland window matching) | `desktop-file-validate` (also in the release workflow) |
| **AppStream** (software centres, Flathub) | `io.github.Shao_RMSGAL.AMCSET.metainfo.xml`: summary, description, license (GPL-3.0-or-later AND the SRIM tables' non-commercial terms), launchable, provides, content rating | `appstreamcli validate`: passed (pedantic: one note). Screenshots are not in it yet (Flathub requires them) |
| **Icon theme specification** | hicolor PNG 16-512 and scalable SVG, named after the application ID | install listing |
| **Manual pages** | `amcset(1)`, `amcset-gui(1)` (version and date from the build) | `groff -ww` without warnings; `man -l` |
| **Packages** | `.deb` (dependencies from `dpkg-shlibdeps`), `.rpm` (automatic requires, `%{dist}` release), AppImage, CLI tarball (C++ runtime linked in, glibc 2.35) | `.deb` built, installed into a folder and run in WSL; the others are built by the release workflow |

## Not done

- Code signing (SignPath Foundation, not set up; 1.0 is unsigned), and the
  App Certification Kit run.
- Screenshots in the AppStream data (and hosting them) before a Flathub
  submission; the Flatpak manifest is unbuilt.
- macOS: no packages (Apple's compiler is Clang, which AMCSET does not
  support yet: AGENTS.md §3).

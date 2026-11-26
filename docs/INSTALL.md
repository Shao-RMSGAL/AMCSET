# Installing AMCSET

AMCSET 1.0 comes as a graphical program (`amcset-gui`) and a command-line
program (`amcset`), with the SRIM-2013 stopping tables (free for
non-commercial use, `data/srim2013/README.md`). Packages are attached to each
release on GitHub (https://github.com/Shao-RMSGAL/AMCSET/releases), with a
`SHA256SUMS.txt` and a build attestation for every file ("Checking a
download" below). To build from source, see `docs/DEVELOPMENT.md`.

| System | Package | Contents |
|---|---|---|
| Windows 10/11 x64 | `AMCSET-<version>-windows-x64.msi` (recommended) | GUI and CLI in Program Files, Start menu entry, `amcset` on the PATH; removable from Settings > Apps |
| | `AMCSET-<version>-windows-x64.zip` | The same files, portable: unzip anywhere, no installation |
| | `AMCSET-<version>-windows-x64.msix` | The same, as an MSIX package (`amcset` as an app alias); only in signed releases (1.0 is not signed) |
| Linux x86-64, any recent distribution | `AMCSET-<version>-x86_64.AppImage` | GUI and CLI in one file |
| Debian 13, Ubuntu 25.04 and later | `amcset_<version>_debian13_amd64.deb` | GUI and CLI, using the system's Qt |
| Fedora 42 and later | `amcset-<version>-1.fc42.x86_64.rpm` | GUI and CLI, using the system's Qt |
| Linux x86-64 (clusters, servers) | `amcset-cli-<version>-linux-x86_64.tar.gz` | The CLI only (no Qt needed) |

## Windows

**MSI.** Double-click it, or install silently from an administrator console:

```powershell
msiexec /i AMCSET-1.0.0-windows-x64.msi /qn
```

It installs to `C:\Program Files\AMCSET` (`bin\amcset-gui.exe`,
`bin\amcset.exe`, `share\amcset\srim2013`, `share\doc\amcset`), adds a Start
menu entry and puts `bin` on the system PATH (open a new console to use
`amcset`). Uninstall from Settings > Apps, or `msiexec /x <file>.msi /qn`.
A newer MSI replaces an older version.

**ZIP.** Unzip, then run `bin\amcset-gui.exe`; `bin\amcset.exe` runs from a
console in that folder. Nothing is written outside the folder except the
GUI's settings and log (`%LOCALAPPDATA%\AMCSET`).

### Running the unsigned Windows release

**AMCSET 1.0 for Windows is not signed.** (Signing through the SignPath
Foundation is planned, `docs/RELEASING.md` §3; a release's notes say whether
its files are signed. A signed file names its publisher under right-click >
Properties > Digital Signatures.) Windows judges a downloaded program by its
signature, so it will question or block these files. That is Windows doing
its job, not a sign that the files are damaged; what tells you the download
is the one this project built is its checksum and build attestation: check
them first ("Checking a download" below).

1. **"Windows protected your PC"** (Microsoft Defender SmartScreen, when the
   MSI or a program from the ZIP is started): choose *More info*, then *Run
   anyway*. To avoid the question, remove the download mark before
   installing or unzipping: right-click the downloaded `.msi` or `.zip` >
   *Properties* > tick *Unblock* > *OK*, or in PowerShell

   ```powershell
   Unblock-File AMCSET-1.0.0-windows-x64.zip
   ```

   Unblock the ZIP **before** unzipping it: files unzipped from a marked
   archive carry the mark too.
2. **Smart App Control** (Windows 11) blocks the programs outright. It lets
   a program run only if Microsoft's cloud service knows it as safe or it
   carries a trusted signature, and Microsoft's FAQ says: "There is
   currently no way to bypass Smart App Control protection for individual
   apps." Its state is under *Settings > Privacy & security > Windows
   Security > App & browser control > Smart App Control settings*. If it is
   **On**, the choices are:
   - turn Smart App Control off there. It then no longer protects against
     any untrusted program, so this is a decision for the whole PC, not for
     AMCSET alone. Microsoft's FAQ says that recent Windows updates allow it
     to be turned on again without reinstalling Windows; on older builds
     switching it off was final until Windows was reset;
   - or leave it on and run the Linux packages in the Windows Subsystem for
     Linux (WSL 2 runs them in its own Linux system; the GUI opens through
     WSLg): "Linux" below;
   - or wait for a signed release.

   The same applies to an AMCSET you build yourself: your own executables
   are unsigned too.
3. **A work or university PC** may forbid unsigned programs by policy (App
   Control for Business, AppLocker). Only its administrators can allow
   AMCSET there; the ZIP needs no administrator rights but is subject to the
   same policy.

Microsoft's descriptions: "What is Smart App Control?"
(https://support.microsoft.com/en-us/topic/what-is-smart-app-control-285ea03d-fa88-4d56-882e-6698afdb7003),
"Smart App Control overview"
(https://learn.microsoft.com/en-us/windows/apps/develop/smart-app-control/overview),
"Microsoft Defender SmartScreen overview"
(https://learn.microsoft.com/en-us/windows/security/operating-system-security/virus-and-threat-protection/microsoft-defender-smartscreen/).

## Linux

**AppImage.** `chmod +x AMCSET-1.0.0-x86_64.AppImage` and run it: the GUI
starts. The CLI is inside: `./AMCSET-1.0.0-x86_64.AppImage amcset run
scenario.json`. It needs FUSE 2 (`libfuse2`) or, without it, `--appimage-extract-and-run`.

**Debian, Ubuntu.** `sudo apt install ./amcset_1.0.0_debian13_amd64.deb`
(it pulls in the Qt libraries it needs). Remove with `sudo apt remove
amcset`.

**Fedora.** `sudo dnf install ./amcset-1.0.0-1.fc42.x86_64.rpm`. Remove with
`sudo dnf remove amcset`.

**CLI only.** Unpack `amcset-cli-1.0.0-linux-x86_64.tar.gz` anywhere and run
`bin/amcset` (the stopping tables are found in `share/amcset` next to
`bin`). It needs only the C and C++ runtime of the system (glibc 2.35 or
later).

The packages install `amcset` and `amcset-gui` to `/usr/bin`, the stopping
tables to `/usr/share/amcset`, the documents to `/usr/share/doc/amcset`, a
desktop entry, AppStream metadata, icons and the manual pages `amcset(1)`
and `amcset-gui(1)`.

## Checking a download

Each release lists the SHA-256 of every file in `SHA256SUMS.txt`:

```powershell
Get-FileHash AMCSET-1.0.0-windows-x64.msi -Algorithm SHA256
```

```bash
sha256sum -c SHA256SUMS.txt --ignore-missing
```

Every file also has a build attestation (GitHub's signed record of the
workflow run that built it from the repository); with the GitHub CLI:

```bash
gh attestation verify AMCSET-1.0.0-x86_64.AppImage --repo Shao-RMSGAL/AMCSET
```

## Where AMCSET looks for its data

The stopping tables are found in this order: the folder named by the
environment variable `AMCSET_DATA_DIR`; `data/` next to the executable or one
or two levels above it; `share/amcset/` next to the executable's folder (the
installed layout); the source tree it was built from.

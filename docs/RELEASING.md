# Releasing AMCSET

How a release is built, signed and published. What the packages contain and
how users install them: `docs/INSTALL.md`. Which platform requirements they
meet: `docs/PLATFORM.md`.

## 1. Packages

| Package | Built by | Where |
|---|---|---|
| `AMCSET-<v>-windows-x64.msi` | CPack + WiX Toolset (`packaging/CMakeLists.txt`, `packaging/windows/wix_patch.xml`) | release workflow, or `cpack -G WIX` locally |
| `AMCSET-<v>-windows-x64.zip` | CPack | release workflow, or `cpack -G ZIP` |
| `AMCSET-<v>-windows-x64.msix` | `tools/windows/make_msix.ps1` (Windows SDK `makeappx`, `packaging/windows/AppxManifest.xml.in`) | release workflow when signing is set up, or locally |
| `AMCSET-<v>-x86_64.AppImage` | `tools/linux/make_appimage.sh` (linuxdeploy) | release workflow (Ubuntu 22.04, Qt 6.8), or a Linux machine |
| `amcset-cli-<v>-linux-x86_64.tar.gz` | CPack, CLI only, `-DAMCSET_STATIC_LIBSTDCXX=ON` | release workflow |
| `amcset_<v>_debian13_amd64.deb` | CPack in a Debian 13 container, system Qt | release workflow |
| `amcset-<v>-1.fc42.x86_64.rpm` | CPack in a Fedora 42 container, system Qt | release workflow |
| `SHA256SUMS.txt`, build attestations | the workflow's release job | release workflow |

Everything comes from `cmake --install`: `bin/` (programs; on Windows also the
Qt and C++ runtime DLLs), `share/amcset/srim2013/`, `share/doc/amcset/`, and on
Linux the desktop entry, AppStream metadata, icons and manual pages
(`packaging/linux/`). The logo and every icon file are drawn by
`scripts/make_icons.py` into `resources/icons/` (committed; rerun it only to
change the logo).

## 2. Cutting a release

1. Version: `project(amcset VERSION x.y.z)` in `CMakeLists.txt` (the programs,
   the Windows version resources, the MSIX, the man pages and the AppStream
   data take it from there) and `version` in `conanfile.py`. Rename the
   CHANGELOG section `## x.y.z (unreleased)` to the release date. Update the
   `<release>` description in
   `packaging/linux/io.github.Shao_RMSGAL.AMCSET.metainfo.xml.in`.
2. Check locally: `tools\windows\build_and_test.bat`, the TRIM comparison if
   the physics changed (`docs/DEVELOPMENT.md` §6), `python3
   scripts/check_headers.py`.
3. Commit, then tag and push: `git tag -a vx.y.z -m "AMCSET x.y.z"`, `git push
   origin vx.y.z`.
4. The workflow `.github/workflows/release.yml` builds and tests on Windows,
   Ubuntu 22.04, Debian 13 and Fedora 42 (about an hour; GitHub-hosted
   runners, free for public repositories), signs the Windows files if signing
   is set up (§3), writes `SHA256SUMS.txt` and build attestations (only in
   a public repository: a private one needs GitHub Enterprise Cloud), and
   creates a **draft** release with every file.
5. Download a few packages from the draft and try them (MSI install and
   uninstall, the AppImage, the CLI tarball on a cluster). Write the release
   notes (the CHANGELOG section), say whether the Windows files are signed,
   and publish the draft.

The workflow can also be started by hand (Actions > Release > Run workflow):
it then builds every package as workflow artifacts, without a release.

**Locally on the workstation** (Windows): after `build_and_test.bat`,

```bat
cd build\Release
cpack -G ZIP -B ..\packages
powershell -ExecutionPolicy Bypass -File ..\..\tools\windows\make_msix.ps1
```

(`cpack -G WIX` needs the WiX Toolset: `dotnet tool install --global wix`
and `wix extension add --global WixToolset.UI.wixext`.) The CPack and CMake
of the Conan build are in `build\Release\generators\conanbuild.bat`'s PATH.
While the GUI is running, `amcset-gui.exe` cannot be replaced (LNK1104) and the
package takes the old one. On Linux: `cpack -G "TGZ;DEB;RPM"` in a build
folder, `tools/linux/make_appimage.sh` for the AppImage.

## 3. Code signing (Windows)

**Why.** Windows judges downloaded programs by their signature. Microsoft
Defender SmartScreen warns about unsigned files and files from publishers
without reputation; **Smart App Control** (on new Windows 11 installations
it starts in evaluation mode and may turn itself on) blocks unsigned programs
outright. Smart App Control accepts
programs signed with a certificate from a certificate authority in the
Microsoft Trusted Root Program (RSA keys; it does not accept ECC
certificates): Microsoft Learn, "Code signing for Smart App Control"
(https://learn.microsoft.com/en-us/windows/apps/develop/smart-app-control/code-signing-for-smart-app-control)
and "Smart App Control overview"
(https://learn.microsoft.com/en-us/windows/apps/develop/smart-app-control/overview).
Every executable and DLL AMCSET ships must therefore be signed (Microsoft's
own runtime DLLs already are), as well as the MSI and the MSIX.

**Decision (maintainer, 2026-11-25): the SignPath Foundation, later.** 1.0
is published unsigned; `docs/INSTALL.md` "Running the unsigned Windows
release" tells users what Windows will do and what they can do about it, and
the draft release's notes say that the Windows files are not signed. §3.2
lists what SignPath needs when the time comes.

Getting a certificate needs an identity check of the publisher and,
usually, payment: that is the maintainer's (or the university's) decision
and action. The options, from Microsoft's "Code signing options for Windows
app developers"
(https://learn.microsoft.com/en-us/windows/apps/package-and-deploy/code-signing-options)
and the providers' pages:

| Option | Publisher shown | Cost and conditions | Fits |
|---|---|---|---|
| **Azure Artifact Signing** (formerly Trusted Signing), Microsoft's signing service | the validated individual or organization | an Azure subscription, identity validation; Basic tier $9.99 per month; individuals in the USA and Canada only at the time of writing, organizations in more countries (https://learn.microsoft.com/en-us/azure/artifact-signing/quickstart) | the release workflow is already written for it (§3.1); certificates are trusted by Smart App Control; keys never leave Microsoft's HSMs |
| **SignPath Foundation** (free for open source) | "SignPath Foundation" | an OSI-approved license (GPL-3.0 is), a public repository, builds on GitHub Actions through SignPath's GitHub integration, multi-factor authentication for the team, an attribution on the project page (https://signpath.org/, https://about.signpath.io/product/open-source) | **the maintainer's choice** (not set up yet, §3.2); needs SignPath's own GitHub Action in the workflow instead of §3.1 |
| **Texas A&M University's code-signing certificate**, if the university has one for its software | Texas A&M University | ask TAMU IT | the publisher name users would expect; sign with `tools\windows\sign.ps1 -Thumbprint` on a machine that holds it |
| An OV code-signing certificate from a commercial CA | the validated organization | a yearly fee; the key on a hardware token or the CA's cloud HSM | `sign.ps1 -Thumbprint` or `-Pfx` locally |

Bundled data note: the SRIM-2013 tables are free for non-commercial use;
that does not affect signing, but the release notes should say so (the
packages' metadata do).

### 3.1 Setting up Azure Artifact Signing for the release workflow

1. In the Azure portal (a subscription is needed): create an **Artifact
   Signing account** in a region (e.g. East US: endpoint
   `https://eus.codesigning.azure.net`).
2. **Identity validation** for the publisher: an individual (USA, Canada) or
   an organization (its legal documents; for the university this goes
   through whoever may act for it). The validated name becomes the
   certificate's subject: what Windows shows as the publisher.
3. Create a **certificate profile** of type Public Trust.
4. Create an **app registration** (Microsoft Entra ID) for GitHub with a
   **federated credential** for this repository (entity: tag, or the
   environment you use), and give it the role **Artifact Signing
   Certificate Profile Signer** on the account. OpenID Connect: no secret is
   stored at GitHub.
5. In the GitHub repository settings: secrets `AZURE_CLIENT_ID`,
   `AZURE_TENANT_ID`, `AZURE_SUBSCRIPTION_ID`; variables
   `ARTIFACT_SIGNING_ENDPOINT`, `ARTIFACT_SIGNING_ACCOUNT`,
   `ARTIFACT_SIGNING_PROFILE`, and `MSIX_PUBLISHER` (the certificate's
   subject exactly as on the profile, e.g. `CN=Texas A&M University, O=...,
   L=College Station, S=Texas, C=US`; the MSIX manifest's publisher must
   equal it).
6. The next tag builds signed files. Check one: right-click > Properties >
   Digital Signatures, or `signtool verify /pa /v <file>`.

The workflow signs with `tools/windows/sign.ps1`: Windows SDK `signtool`
with Microsoft's signing client (NuGet package
`Microsoft.ArtifactSigning.Client`, `Azure.CodeSigning.Dlib.dll`; needs the
.NET 8 runtime), SHA-256 digests and an RFC 3161 time stamp from
`http://timestamp.acs.microsoft.com` (Artifact Signing's certificates last
three days; the time stamp keeps signatures valid after that): Microsoft
Learn, "Set up signing integrations to use Artifact Signing"
(https://learn.microsoft.com/en-us/azure/artifact-signing/how-to-signing-integrations).
It signs `amcset.exe`, `amcset-gui.exe` and the Qt DLLs in the build folder
before CPack packs them, then the MSI and the MSIX; files that already have
a valid signature are left alone.

**Locally**, with the same account (after `az login` as an identity that has
the signer role), or with a certificate in the Windows store:

```powershell
powershell -File tools\windows\sign.ps1 -Metadata metadata.json -Dlib <client>\bin\x64\Azure.CodeSigning.Dlib.dll build\packages
powershell -File tools\windows\sign.ps1 -Thumbprint <thumbprint> -TimestampUrl <CA's RFC 3161 URL> build\packages
```

SmartScreen reputation builds up over the first downloads of files signed
with a certificate; until then a warning may still appear for signed files.

### 3.2 SignPath Foundation (chosen; to be set up)

Nothing of this exists yet. What it takes, from SignPath's pages
(https://signpath.org/, https://about.signpath.io/product/open-source);
read their current terms before applying, they decide:

1. The maintainer applies for the project at signpath.org (the license, the
   public repository, who the team members are, multi-factor authentication
   on their accounts). Creating that account and accepting SignPath's terms
   are the maintainer's actions.
2. SignPath's conditions for the project page and the releases (an
   attribution naming SignPath, a code-signing policy saying who may approve
   a release) go into `README.md` and this document.
3. In `release.yml`, the Windows job uploads the unsigned build as an
   artifact and SignPath's GitHub Action submits a signing request for it
   (SignPath signs what a GitHub-hosted run built from this repository, not
   files from a workstation); the signed files come back as an artifact and
   are packed and attached as now. The Azure steps (§3.1) run only when
   `ARTIFACT_SIGNING_ENDPOINT` is set, so they can stay unused.
4. The MSIX needs `MSIX_PUBLISHER` to equal the subject of SignPath's
   certificate; until then it is not built.
5. Check a signed file with `signtool verify /pa /v <file>`; the release
   notes of that release then say that the Windows files are signed.

## 4. Linux: checks instead of signatures

Linux packages from a GitHub release are checked by their SHA-256
(`SHA256SUMS.txt`) and by GitHub's build attestations (`gh attestation verify
<file> --repo Shao-RMSGAL/AMCSET`: the file was built by this repository's
release workflow from the tagged commit). Optionally, sign `SHA256SUMS.txt`
with the maintainer's GPG key after downloading it from the draft (`gpg
--armor --detach-sign SHA256SUMS.txt`) and attach the `.asc`. Distribution
repositories (Debian, Fedora, Flathub) sign their own packages.

## 5. Catalogues (after the first release)

- **winget** (Windows Package Manager): `wingetcreate new
  <URL of the MSI in the release>` writes the manifests from the installer
  and opens a pull request to `microsoft/winget-pkgs`; for later versions
  `wingetcreate update`. A signed MSI is strongly advisable (winget's
  validation scans the installer).
- **Microsoft Store**: the signed MSIX (or the MSI) through Partner Center;
  needs a developer account. The Store re-signs MSIX packages itself.
- **Flathub**: `packaging/flatpak/io.github.Shao_RMSGAL.AMCSET.yml` is a
  starting manifest (KDE runtime 6.8; abseil and Boost built in). It has not
  been built yet: fill in the source checksums, build with `flatpak-builder`,
  then follow Flathub's submission process. Flathub must be told that the
  bundled stopping tables are non-commercial (`project_license` in the
  AppStream data says so).

## 6. Before publishing

- [ ] The CHANGELOG section, the version, the AppStream release text.
- [ ] The workflow passed on every platform (tests included).
- [ ] An MSI installed, `amcset version` in a new console, the GUI from the
      Start menu, uninstalled from Settings > Apps.
- [ ] The AppImage and the `.deb` or `.rpm` start the GUI; the CLI tarball runs.
- [ ] Windows files signed, or (1.0) the release notes say they are not and
      point to docs/INSTALL.md "Running the unsigned Windows release" (the
      workflow writes that sentence when signing is not configured).
- [ ] `SHA256SUMS.txt` attached; the release notes link `docs/INSTALL.md`.

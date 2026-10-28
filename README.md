# AMCSET

AMCSET is a C++ application for simulating ion and electron bombardment of
layered materials. Its desktop interface lets you configure a target and
bombardment, run a simulation, and inspect the resulting particle behavior and
summary data. It is intended as a tool for exploring material response to
particle bombardment, including energy loss, displacement, and particle range.

The application is under active development. The graphical application is the
primary interface; the command-line mode is not yet implemented.

## How it works

The project is split into two main parts:

- `common/` contains the simulation and materials code shared by the
  application. It models particle interactions and uses bundled isotope and
  scattering data.
- `gui/` is a Qt 6 desktop application for setting up layered targets and
  simulations, visualizing results, and configuring output.

The simulation can collect results such as displacement, electronic and
nuclear energy-loss, total energy-loss, and range histograms. The GUI allows
you to choose output locations; its default output directory is `output/`.
The project uses C++20, CMake, and Conan to manage its C++ dependencies.

## Build and run on Linux

Ubuntu 24.04 is the environment used by the project's Docker image. Building
Qt and other dependencies from source can take a long time, especially on the
first build.

Install the build tools and Conan 2.x. You will need a C++20-capable compiler
(GCC 14 is used in the Docker image), CMake 3.23 or newer, Ninja, Python 3,
pipx, and the system development libraries required by Qt 6 (including the
X11/Wayland and OpenGL/Vulkan development libraries). On Ubuntu, the base
tools can be installed with:

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build python3 python3-pip python3-venv pipx
pipx install "conan>=2,<3"
export PATH="$HOME/.local/bin:$PATH"
```

If your distribution does not provide a recent enough CMake, Ninja, compiler,
or Qt system dependency, use the Docker workflow below. Then, from the
repository root:

```bash
conan profile detect --force
./qbuild.sh -c -b
```

`-c` installs the dependencies through Conan and `-b` configures and builds
the Release application. The first dependency install may compile packages
from source and take substantial time. The detected Conan profile must match
your installed compiler; the build requests the GNU C++20 and C17 standards.

Once built, use the helper script to run tests or launch the GUI:

```bash
./qbuild.sh -t  # Run the tests
./qbuild.sh -r  # Launch the GUI
```

You can also perform all three steps in sequence with
`./qbuild.sh -c -b -t -r`. The build and test files are placed under
`build/Release/`.

## Build and run in Docker

The Docker image provides an Ubuntu 24.04 development environment with pinned
compiler and build-tool versions. The included `run-dev.sh` helper is designed
for Linux hosts: it builds the image if needed, mounts this repository into
the container, persists the Conan cache in a Docker volume, and forwards an
available Wayland or X11 display (and `/dev/dri` when available).

Requirements:

- Docker Engine and the Docker CLI
- A Linux desktop session for displaying the GUI
- For VS Code integration, the Dev Containers extension

From the repository root, start the development container:

```bash
./run-dev.sh
```

The first image build and Conan dependency setup can take a while. When the
script finishes, attach VS Code to `amcset-dev-container` using **Dev
Containers: Attach to Running Container**. In a terminal in the container,
build, test, and run as usual:

```bash
./qbuild.sh -b
./qbuild.sh -t
./qbuild.sh -r
```

`run-dev.sh --rebuild` forces the Docker image to be rebuilt. The workspace is
mounted from the host, so source changes and generated build files remain in
the repository. Docker keeps downloaded Conan packages in a named volume so
they can be reused on later runs. If the host has no usable display, the
container can still be used to build and test, but the GUI will not appear.

## Windows

The PowerShell helper builds the Release configuration using Conan, CMake, and
the compiler toolchain configured in your Conan profile. Install those tools
first (including Python and a supported C++20 compiler), then open PowerShell
in the repository root:

```powershell
.\qbuild.ps1 -c  # Install Conan dependencies
.\qbuild.ps1 -b  # Configure and build
.\qbuild.ps1 -t  # Run tests
.\qbuild.ps1 -r  # Launch the GUI
```

The dependency build may take a long time on the first run. The script supports
the `-c`, `-b`, `-t`, and `-r` options; install prerequisites separately.

## License

AMCSET is distributed under the GNU General Public License, version 3 or (at
your option) any later version. See [LICENSE](LICENSE).

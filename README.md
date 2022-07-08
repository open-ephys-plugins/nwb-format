# NWB Format

![header-image](Resources/header-image.png)

A Record Engine plugin for the Open Ephys GUI based on the [Neurodata Without Borders 2.X specification](https://nwb-schema.readthedocs.io/en/latest/format_release_notes.html).

## Installation

This plugin can be added via the Open Ephys GUI's built-in Plugin Installer. First, make sure there are no Record Nodes in your signal chain. Then, press **ctrl-P** or **⌘P** to open the Plugin Installer, browse to the "NWB Format" plugin, and click the "Install" button. The next time you add a Record Node to the signal chain, `NWB2` will be available in the data format drop-down menu.

## Usage

The specifications of NWB files written by the Open Ephys GUI are available [here](https://open-ephys.github.io/gui-docs/User-Manual/Recording-data/NWB-format.html).

## Building from source

First, follow the instructions on [this page](https://open-ephys.github.io/gui-docs/Developer-Guide/Compiling-the-GUI.html) to build the Open Ephys GUI.

This plugin depends on the [OpenEphysHDF5Lib](https://github.com/open-ephys-plugins/OpenEphysHDF5Lib) common library. Make sure you build and install that first before you proceed with building NWBFormat plugin.

**Important:** This plugin is intended for use with the pre-release core application, version 0.6.0. The GUI should be compiled from the [`development-juce6`](https://github.com/open-ephys/plugin-gui/tree/development-juce6) branch, rather than the `master` branch.

Then, clone this repository into a directory at the same level as the `plugin-GUI`, e.g.:
 
```
Code
├── plugin-GUI
│   ├── Build
│   ├── Source
│   └── ...
├── OEPlugins
│   └── nwb-format
│       ├── Build
│       ├── Source
│       └── ...
```

### Windows

**Requirements:** [Visual Studio](https://visualstudio.microsoft.com/) and [CMake](https://cmake.org/install/)

From the `Build` directory, enter:

```bash
cmake -G "Visual Studio 17 2022" -A x64 ..
```

Next, launch Visual Studio and open the `OE_PLUGIN_nwb-format.sln` file that was just created. Select the appropriate configuration (Debug/Release) and build the solution.

Selecting the `INSTALL` project and manually building it will copy the `.dll` and any other required files into the GUI's `plugins` directory. The next time you launch the GUI from Visual Studio, `NWB2` should appear as an data format option in the Record Node.


### Linux

**Requirements:** [CMake](https://cmake.org/install/)

From the `Build` directory, enter:

```bash
cmake -G "Unix Makefiles" ..
cd Debug
make -j
make install
```

This will build the plugin and copy the `.so` file into the GUI's `plugins` directory. The next time you launch the compiled version of the GUI, `NWB2` should appear as an data format option in the Record Node.


### macOS

**Requirements:** [Xcode](https://developer.apple.com/xcode/) and [CMake](https://cmake.org/install/)

From the `Build` directory, enter:

```bash
cmake -G "Xcode" ..
```

Next, launch Xcode and open the `nwb-format.xcodeproj` file that now lives in the “Build” directory.

Running the `ALL_BUILD` scheme will compile the plugin; running the `INSTALL` scheme will install the `.bundle` file to `/Users/<username>/Library/Application Support/open-ephys/plugins-api`. `NWB2` should now appear as an data format option in the Record Node.


### Attribution

This plugin, along with the OpenEphysHDF5Lib, were developed collaboratively by Aaron Cuevas Lopez, Pavel Kulik, and Josh Siegle. It is based on the [NWB Format Specification](https://nwb-schema.readthedocs.io/en/latest/format_release_notes.html).
 

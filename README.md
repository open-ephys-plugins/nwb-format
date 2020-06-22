# NWBFormat
A format to record to files based on the Neurodata Without Borders 1.0 specification

## Installation
### Installing the OpenEphysHDF5Lib Common Library
This plugin depends on the [OpenEphysHDF5Lib](https://github.com/open-ephys-plugins/OpenEphysHDF5Lib) common library. Make sure you build and install that first before you proceed with building NWBFormat plugin.

### Building the plugins
Building the plugins requires [CMake](https://cmake.org/). Detailed instructions on how to build open ephys plugins with CMake can be found in [our wiki](https://open-ephys.atlassian.net/wiki/spaces/OEW/pages/1259110401/Plugin+CMake+Builds).
We highly recommend building all three projects simultaneously using the configuration provided in the top-level Build folder, instead of the individual plugins. Building the plugins individually will need manual tweaking for them to find the OpenEphysHDF5 common library.
 
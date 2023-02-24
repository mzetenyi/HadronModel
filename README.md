# HadronModel
Codes for calculations in an Effective Lagrangian model for hadron reactions.

The project is organized in two libraries and several executables. The libraries are:

- **`FeynTools`**: A collection of general tools useful for calculating feynman digrams in relativistic field theory. This includes representation of half-integer numbers, Lorenty fourvectors and fourtensors, Lorentz transformations, Dirac spinors and gamma matrices, etc. together with some utilities to e.g. configure global parameters, handle measuring units.
- **`ELModel`**: The implementation of an Effective Lagrangian model used for calculating di-electron and pion-pair production in pion-nucleon collisions in the resonance region.

Some of the executables are:

- **`piN_Ndilep`**: A code for calculating pion induced di-electron production, $\pi + N \to N + e^+ + e^-$. There are various outpus for observables like differential cross section $d\sigma/d M_{e^+e^-}$, and spin density matrix elements of the virtual photon decaying to the lepton pair.
- **`piN_Npipi`**: A code for calculating contributions of a $\rho$ meson to pion pair production in pion-nucleon collisions, $\pi + N \to N +\rho \to N + \pi + \pi$.
- **`piN_Npipi_generator`**: Event generator for the process $\pi + N \to N +\rho \to N + \pi + \pi$.
- **`NN_NDelta`**: Code for calculating the differential cross section $d\sigma/dm_\Delta$ of the reaction $N+N \to N+\Delta$.

## How to download and compile

### Download
Navigate to the [main page of the project](https://github.com/mzetenyi/HadronModel), click on the green "Code" button and choose **"Download zip"**. Save the zip file to a suitable directory, unzip it, and you will get a directory structure like this:

![Directory structure](https://user-images.githubusercontent.com/43382422/220627571-33fd9b92-198a-4633-8c0c-6331880fbb33.jpg)

### Compile

#### **prerequisites**
You will need `cmake` version 3.0.0 or higher, and a c++17 compliant compiler. For the developement, `cmake` version 3.16.3 and `gcc` version 9.4.0 was used in a linux environment. 
#### **steps of compilaton**
Navigate to the main directory of the project, create a build directory, and navigate to it:
```
> mkdir build
> cd build
```
Generate the build system:
```
> cmake ..
```
Generate a specific executable:
```
> cmake --build . --target <target_name>
```
where `<target_name>` is the name of one of the executables listed above, e.g. `piN_Npipi_generator`. You can omit the `--target` specification to build all targets:
```
> cmake --build .
```
During compilation some warnings are issued which you can ignore.

In the build directory, the folder structure of the project is repeated with source files replaced by the corresponding executables. E.g. the executable for the $\pi\pi$ event generator can be found here:
```
HadronModel-main/build/exec/pi-N/piN_Npipi_generator
```
The dirctory structure of the project looks like the following after compilation:
![Directiries2](https://user-images.githubusercontent.com/43382422/221149994-5b1ddf68-125d-41d6-9bbd-c9dcff99e071.jpg)

## Executables

Executables are run from the command line and write the results to standard output that the user should redirect to a file. Parameters are specified via various command line options.

There are two types of options:
1. flags, e.g. for chosing the output type or swiching reaction channels on and off
2. parameters with values given in the `<key>=<value>` format (spaces around the = sign are NOT allowed!).

The special option `load[<path_to_options_file>]` loads the given file that should contain one option specification per line, an example is [the file specifying the model parameters](https://github.com/mzetenyi/HadronModel/blob/main/lib/ELModel/model_params). Observe the comment syntax `//`, and the possibility for grouping the options using the syntax
```
<gruop> {
  <key1> = <value1>
  <key2> = <value2>
}
```
Options in files can be overridden on the command line, for grouped parameters use the syntax `<group>.<key1>=<value1>`. Command line options are read in from left to right, options read in later always override previously read-in values of the same option.

### piN_Npipi_generator

This executable is an event generator for contributions of a $\rho$ meson to pion pair production in pion-nucleon collisions, i.e. the process $\pi^- + p \to n +\rho \to n + \pi^+ + \pi^-$.

### Options
#### **Options with values:**

**srt:** The CM energy, $\sqrt{s}$ given in GeV.

**Ngen:** The number of generated events.

**res_rhophase:** An extra phase factor for all baryon resonance contributions with respect to the Born term. Given in degrees.

#### **Flags:**

**flags for channels:** The following options can be used to switch on the corresponding reaction channel: `Born`, `N1440`, `N1520` , `N1535`, `N1650`, `N1675`, `N1680`, `N1710`, `D1600`, `D1620`.

A typical command to run the code, issued in the main directory of the project:
```
./build/exec/pi-N/piN_Npipi_generator load[lib/ELModel/model_params] srt=1.49 res_rhophase=90 Born N1520 Ngen=1000 > my_output
```
This will generate 1000 events of pion pair production at CM energy $\sqrt{s}=1.49\ \text{GeV}$, including the Born and $N(1520)$ channels with an extra relative phase of $90^\circ$, and write the results in the file `my_output`.

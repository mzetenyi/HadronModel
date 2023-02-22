# HadronModel
Codes for calculations in an Effective Lagrangian model for hadron reactions.

The project is organized in two libraries and several executables. The libraries are:

- **FeynTools**: A collection of general tools useful for calculating feynman digrams in relativistic field theory. This includes representation of half-integer numbers, Lorenty fourvectors and fourtensors, Lorentz transformations, Dirac spinors and gamma matrices, etc. together with some utilities to e.g. configure global parameters, handle measuring units.
- **ELModel**: The implementation of an Effective Lagrangian model used for calculating di-electron and pion-pair production in pion-nucleon collisions in the resonance region.

Some of the executables are:

- **piN_Ndilep**: A code for calculating pion induced di-electron production, $\pi + N \to N + e^+ + e^-$. There are various outpus for observables like differential cross section $d\sigma/d M_{e^+e^-}$, and spin density matrix elements of the virtual photon decaying to the lepton pair.
- **piN_Npipi**: A code for calculating contributions of a $\rho$ meson to pion pair production in pion-nucleon collisions, $\pi + N \to N +\rho \to N + \pi + \pi$.
- **piN_Npipi_generator**: Event generator for the process $\pi + N \to N +\rho \to N + \pi + \pi$.
- **NN_NDelta**: Code for calculating the differential cross section $d\sigma/dm_\Delta$ of the reaction $N+N \to N+\Delta$.

## How to download and compile

### Download
Navigate to the [main page of the project](https://github.com/mzetenyi/HadronModel), click on the green "Code" button and choose **"Download zip"**. Save the zip file to a suitable directory, unzip it, and you will get a directiry structure like this:

![Directories](https://user-images.githubusercontent.com/43382422/220627571-33fd9b92-198a-4633-8c0c-6331880fbb33.jpg)

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

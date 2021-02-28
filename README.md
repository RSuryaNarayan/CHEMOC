# What is CHEMOC
CHEMOC is an acronym for CHEMically reacting Method Of Characteristics. This is an attempt to implement solutions to governing PDEs of suspersonic non-equilibrium reacting gas-dynamics using the classical method of characteristics approach. Method of Characteristics (MOC) has been used previously to seek solutions to supersonic flow fields by marching along characteristics to chart out the flow field along Mach-lines. In CHEMOC we use this approach, however, in-addition perform 1-D finite-rate chemistry to account for chemical reactions. Popular uses of MOC include computing the supersonic nozzle contour, Taylor-Macoll flow around ogival noses and much more. All of these require suitable unit processes that could be used to iterate to produce the flow field via marching along characteristics. CHEMOC is an attempt to provide pre-coded modules for standard 3-characteristic unit-process and make solutions to these flows simple. 
# Dependencies and setup
CHEMOC is written entirely in the Julia programming language. CHEMOC uses Cantera to integrate chemistry into flow physics. Cantera is integrated into Julia using the `PyCall` library. detailed steps to integrate Cantera into Julia is given below. For updates on a native Julia code for Cantera-like packages, follow [Cantera-enhancements](https://github.com/Cantera/enhancements/issues/81). Attempts are being made to benchmark finite-rate chemistry cases with CHEMOC. To start using CHEMOC you first have to integrate Cantera and Julia. 
## Integrating Cantera and Julia
Go to the [Cantera Downloads](https://cantera.org/install/index.html) page and install Cantera in your system using [Python installation](https://cantera.org/install/windows-install.html). Follow the steps in the website as is and test if Cantera works in the Python IDLE. Copy the location of the python.exe. Alternatively you can and download [these](https://drive.google.com/drive/folders/1qF9Au-2RxchwKZ-GAYi-oBdsODER0j4U) and follow theses steps:
* Install python (run the python setup). 
* configure numpy using shell
* run the cantera setup
* configure cantera in python
* copy the location of the python.exe executable
Now open the Julia REPL and install `PyCall`
```
using Pkg;
Pkg.add("PyCall");
```
Once Julia is done updating its registries, update the local python environment that has the Cantera installation using the following commands:
```
using PyCall
ENV["PYTHON"] = "location of the python.exe"
```
The location should look something like `C:\Users\<Username>\AppData\Local\Programs\Python\Python37`,if you are using windows. Now STOP the REPL and restart Julia. You can now run Cantera using: 
```
using PyCall;
pyimport("cantera");
gas1=cantera.Solution("gri30.xml");
```
What we have done is initialized a `gas` object using the `GRI30` combustion mechanism of natural gas with with reburn chemistry. You should get a result like so: 
```
gri30:
temperature    300 K
pressure    101325 Pa
density  0.0818891 kg/m^3
mean mol. weight 2.01588 amu
```
It should also show other properties like `internal energy`, `enthalpy`, `Gibbs function`, `heat capacity c_p`, `heat capacity c_v`.
## Usage
CHEMOC is originally targetted at producing a shock-free nozzle contour. However there exist some generic functions with self-explanotory names to march along specific characteristics. More its usage to be done. 
## Results so far
An initial kernel for the rocket nozzle design problem has been obtained using the `Unit_process_AXIS` and `Unit_process_INTERIOR_Cminus` functions:
![plot!](https://github.com/RSuryaNarayan/CHEMOC/blob/main/Results/Kernel.PNG)\
The following is plot of the Mach number along the axis of the nozzle.
![plot!](https://github.com/RSuryaNarayan/CHEMOC/blob/main/Results/Mach%20number.PNG)

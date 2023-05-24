# PX915 Group Project - Group A

This is a simulation code for a half-cell single particle model (half-SPM model) written in Fortran90 and interfaced with Python. The code calculates the change in concentration in a single, spherical particle for a set simulation time from which a voltage profile can also be obtained.
The relevant input parameters can be set in the Jupyter notebook ‘Input_NetCDF.ipynb’ file, which saves the parameters as a NetCDF file. This will be read by the programme, which runs the simulation and returns a visualisation of the concentration and (if chosen) of the voltage profile.

## Getting started
A Jupyter notebook tutorial is provided to guide the user through the the main features of the simulation package. 
This is most easily accessed by cloning the GitHub repository to the machine that will be used to run the simulation.

```
git clone https://github.com/HubertJN/PX915_Group_Project.git
```

The tutorial notebook can be accessed by moving to the directory containing the cloned repository and opening a new jupyter notebook environment.

```
cd [directory the repository has been cloned to]
jupyter notebook
```
There is the option to either start a new simulation or to continue a simulation from a checkpoint file. In the first case, all the input parameters have to be set and a new input file in NetCDF format will be created. If the simulation is continued from the checkpoint file, all input parameters are kept the same as for the first run except for the simulation length and the number of simultion steps after which a new checkpoint file is created.

## Setting the input parameters
The default parameters are for a positive NCM (Nickel-Cobalt-Manganese) electrode and were taken from [Chen2020](https:/doi.org/10.1149/1945-7111/ab9050).

The notebook provides slider bars and an input cell that allow the user to change the input parameters within suggested ranges. There is also the option so set the parameters in a unrestricted input cell – however the user is advised to use this sensibly, since it otherwise might lead to unphysical behaviour or failure of the code.

Important consideration are:
- time step: if the time steps are set too high, the simulation might take too long to complete
- SOC and applied current: at 0% or 100% state of charge (SOC), the applied current should be set positive (charging) or negative (discharging/use), respectively
- If an anode material is being simulated, it will still be the positive electrode with respect to lithium in this half-cell SMP model and the parameters need to be set accordingly

## Running the simulation
Once the input parameters are set and the 'SPM_input.nc' file has been created, the simulation can be run using the provided Makefile. 
To compile the code the command is:
```
make
```
To run the simulation, the command is:
```
make exe
```
To visualise the output, the command is:
```
make visual
```

## Documention
For more information see our [documentation](./docs.pdf).
If the documentation needs to be re-generated, input the following into the command prompt.
```
make docs
```

## Uncertainty Quantification
There are a variety of option to visualise the uncertainties involved with the simulation. The following commands represent the available uncertainty quantification options.

Perform sensitivity analysis and then display the results.
```
make sensitive
```

Display the results from the sensitivity analysis.
```
make vis_sens
```

Perform uncertainty propagation using random latin hypercube sampling and display the results.
```
make uncertain
```

Display the results from random latin hypercube sampling.
```
make uncer_vis
```

Perform sensitivity analysis and perform random latin hypercube sampling. Then calculate uncertainty from the standard deviations and random latin hypercube sampling. Then display results.
```
make sens_uncer
```

Visualise results of calculated uncertainties from the standard deviations and random latin hypercube sampling.
```
make vis_sens_uncer
```

Calculates the uncertainty from the standard deviations and then, assuming random latin hypercube sampling has been performed, displays both.
```
make sens_uncer_sep
```
## How to run the data fitting
Download the /venv/ directory, which contains all the necessary python libraries which the data fitting notebook, 'DataFit_int.ipynb', requires. The /venv/ directory is a virtual environment, and launching Jupyter notebook inside this virtual environment will allow Jupyter to use the required libraries without the user needing to install them. To activate the virtual environment, go to the directory where /venv/ is stored. This directory should also contain the solver used for data fitting, 'datafitPde.f90', and the notebook that runs the datafitting, 'DataFit_int.ipynb'. Once the /venv/ directory and the aforementioned programs are stored in the same location, the virtual environment needs to be launched. To launch the virtual environment, move to the location where it is stored in a terminal and type the following command:
```
source /venv/bin/activate
```
You can tell if the virtual environment is activated if you can see '(venv)' written before your current working directory in the terminal. Jupyter notebook can now be launched in this virtual environment by typing the following in the same terminal:
```
jupyter notebook
```
Once Jupyter is launched, find the notebook 'DataFit_int.ipynb' in the Jupyter file browser and open it. You can then run the cells accordingly to preform the optimisation on the diffusion coefficient.


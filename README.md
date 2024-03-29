# Half-SPM-Solver

This is a simulation code for a half-cell single particle model (half-SPM model) written in Fortran90 and interfaced with Python. The code calculates the change in concentration in a single, spherically symmetric particle for a set simulation time from which a voltage profile can also be obtained.
The relevant input parameters can be set in the Jupyter notebook ‘Input_NetCDF.ipynb’ file, which saves the parameters as a NetCDF file. This will be read by the program, which runs the simulation and returns a visualisation of the concentration and (if chosen) of the voltage profile.

## Getting Started
A Jupyter notebook tutorial is provided to guide the user through the the main features of the simulation package. 
This is most easily accessed by cloning the GitHub repository to the machine that will be used to run the simulation.
<br>
```
git clone https://github.com/HubertJN/PX915_Group_Project.git
```
In a terminal, move into the directory where the repository has been cloned with the command:
<br>
```
cd [directory the repository has been cloned to]
```
If the software is being used for the first time, the virtual enviroment must be set up beforehand. In the same directory, this can be done with the following commands:
<br>
```
make virtual
source venv/bin/activate
make mods
```
These commands create a folder called "venv/" and downloads all the necessary libraries for the Python programs in the software. This removes the necessity for the user to manually install all the dependencies of the software.

Once the virtual environment is set up for the first time, it must always be activated before using the software. To activate the virtual environment, the following command must be executed:
<br>
```
source venv/bin/activate
```
For an indication that the virtual environment is activated, it can be checked that '(venv)' is written before the current working directory in the terminal. 

The tutorial notebook can be accessed by opening a new jupyter notebook environment:
<br>
```
jupyter notebook
```
Once Jupyter is launched, find the notebook 'Tutorial.ipynb' in the Jupyter file browser and open it.There is the option to either start a new simulation or to continue a simulation from a checkpoint file. In the first case, all the input parameters have to be set and a new input file in NetCDF format will be created. If the simulation is continued from the checkpoint file, all input parameters are kept the same as for the first run except for the simulation length, the number of simultion steps and the applied current after which a new checkpoint file is created.

To deactivate the virtual environment after using the software, the simple following command can do this:
<br>
```
deactivate
```

## Setting the Input Parameters
The default parameters are for a positive NCM (Nickel-Cobalt-Manganese) electrode and were taken from [Chen2020](https:/doi.org/10.1149/1945-7111/ab9050).

The notebook provides slider bars and an input cell that allow the user to change the input parameters within suggested ranges. There is also the option so set the parameters in a unrestricted input cell – however the user is advised to use this sensibly, since it otherwise might lead to unphysical behaviour or failure of the code.

Important consideration are:
- time step: if the time steps are set too high, the simulation might take too long to complete
- SOC and applied current: at 0% or 100% state of charge (SOC), the applied current should be set positive (charging) or negative (discharging/use), respectively
- The convention for this solver when initialising the applied current density is opposite to the general convention. In this solver, a positive applied current density corresponds to a charging battery, and a negative applied current density corresponds to a discharging battery.
- If an anode material is being simulated, it will still be the positive electrode with respect to lithium in this half-cell SMP model and the parameters need to be set accordingly

## Running the Simulation
Once the input parameters are set and the 'SPM_input.nc' file has been created, the simulation can be run using the provided Makefile.

### Serial Code
To run the code in serial execute the following commands:
<br>
```
make
```
To run the simulation, the command is:
<br>
```
make exe
```
To visualise the output, the command is:
<br>
```
make visual
```
### Parallel Code
To run in parallel open "Makefile" and set "num_threads" to the desired number of threads to parallelise over. Run same commands as above.

## Uncertainty Quantification
There are a variety of options to visualise the uncertainties involved with the simulation. The following commands represent the available uncertainty quantification options.

Perform sensitivity analysis and then display the results.
<br>
```
make sensitive
```

Display the results from the sensitivity analysis.
<br>
```
make vis_sens
```

Perform uncertainty propagation using random latin hypercube sampling and display the results.
<br>
```
make uncertain
```

Display the results from random latin hypercube sampling.
<br>
```
make uncer_vis
```

Perform sensitivity analysis and perform random latin hypercube sampling. Then calculate uncertainty from the standard deviations and random latin hypercube sampling. Then display results.
<br>
```
make sens_uncer
```

Visualise results of calculated uncertainties from the standard deviations and random latin hypercube sampling.
<br>
```
make vis_sens_uncer
```

Calculates the uncertainty from the standard deviations and then, assuming random latin hypercube sampling has been performed, displays both.
<br>
```
make sens_uncer_sep
```
## How to Run the Data Fitting

To run the data fitting, first move into the datafitting folder:
<br>
```
cd datafitting/
```
In this datafitting folder, compile the data fitting software with the command:
<br>
```
./compile.sh
```
The data fitting fortran file is compiled so that it can be wrapped into the python jupyter notebook.

Jupyter notebook can now be launched in the virtual environment by typing the following in the same terminal:
<br>
```
jupyter notebook
```
Once Jupyter is launched, find the notebook 'DataFit_int.ipynb' in the Jupyter file browser and open it. You can then run all the cells in order to perform the optimisation on the diffusion coefficient.
To see the result in a separate window open the png file 'fitted_diff_coeffs.png'.

## How to obtain the benchmarking results
The half-SPM model has been benchmarked against the open source simulation package PyBaMM by obtaining the relative absolute error and the root mean square error (RMSE) for 5 different simulations:

1. Charging at 5 A using the default input parameters
2. Discharging at 5 A starting at 95% state of charge (SOC) using the default for all other input parameters
3. A galvanostatic intermittent titration technique (GITT) experiment, where the half-cell is rested for 30 minutes (i.e. no current applied), followed by discharge at 1 A for 5 minutes and another rest period for 30 minutes
4. Charging with a diffusion coefficient of $4*10^{-10} \ m^2 s^{-1}$ using the default for all other input parameters
5. Charging with a mean particle radius of $7 \mu m$ using the default for all other input parameters

A convergence study of the RMSE with respect to the time step is also provided.

The results can be obtained by running the following comamnd:

```
make benchmarking
```

## Documention
For more information see our [documentation](./docs.pdf).
If the documentation needs to be re-generated, input the following into the command prompt.
<br>
```
make docs
```

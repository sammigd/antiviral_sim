# antiviral_sim
Simulation study comparing antiviral treatment strategies for a hypothetical influenza season

Script 0
- Creates some parameter inputs using code from this repo: https://github.com/CDCgov/flu-vaccine-model

Script 1
- Sets parameters that do not change across simulation scenarios

Script 2
- Sets scenario specific parameters (e.g. defines ranges considered for transmissibility, antiviral effectiveness, etc)

Script 3
- Function with ordinary differential equations defining state transitions for compartmental model

Script 4
- Calls the prior scripts to create parameters and model
- Runs model for all parameter sets, generates and saves output

To run the model, just run script 4 in its entirety. If you have an R project or working directory set to the repo home folder, everything should run fine. 

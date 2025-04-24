# Influenza antiviral treatment population model
Simulation study comparing antiviral treatment strategies for a hypothetical influenza season

Code is in the model_code folder.

0_create_contact_matrix.R
- Creates some parameter inputs using code from this repo: https://github.com/CDCgov/flu-vaccine-model
- No one should really need to run this. I didn't add the scripts from the other repo since the needed outputs are saved in the inputs folder already.

1_set_fixed_params.R
- Sets parameters that do not change across simulation scenarios

2_set_scenario_params.R
- Sets scenario specific parameters (e.g. defines ranges considered for transmissibility, antiviral effectiveness, etc)

3_create_model_function.R
- Function with ordinary differential equations defining state transitions for compartmental model

4_generate_output.R
- Calls the prior scripts to create parameters and model
- Runs model for all parameter sets, generates and saves output

To run the model, just run script 4 in its entirety. If you have an R project or working directory set to the repo home folder, everything should run fine. 

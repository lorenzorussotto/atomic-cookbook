import numpy as np
import os
from temp_dict import temps_55, temps_147, temps_309

atoms = [55, 147, 309]
percentages = [0.7, 0.76, 0.8, 0.86, 0.9, 0.96, 1]

timestep = 0.005
steady_state_steps = 2000000

steps_per_day = {
    55: 12400000,
    147: 16000000,
    309: 12200000
}

temp_dicts = {
    55: temps_55,
    147: temps_147,
    309: temps_309
}
  
def compute_sim_features(atoms, pctg):
    # Composition features
    gold_atoms = int(np.trunc(atoms * pctg))
    gold_percentage = int(round(pctg * 100, 0))
    pd_percentage = 100 - gold_percentage
    daily_steps = steps_per_day[atoms]

    return(gold_atoms, gold_percentage, pd_percentage, daily_steps)

def momentum_fix(pctg, script):
    if pctg == 0:
            script = script.replace("${MOMENTA_AU}", "") 
            script = script.replace("${MOMENTA_PD}", "fix             fix1 Pd momentum 1 linear 1 1 1 angular # fix Pd angular momenta")

    elif pctg == 1:
            script = script.replace("${MOMENTA_AU}", "fix             fix2 Au momentum 1 linear 1 1 1 angular # fix Au angular momenta") 
            script = script.replace("${MOMENTA_PD}", "")

    else:
            script = script.replace("${MOMENTA_PD}", "fix             fix1 Pd momentum 1 linear 1 1 1 angular # fix Pd angular momenta")
            script = script.replace("${MOMENTA_AU}", "fix             fix2 Au momentum 1 linear 1 1 1 angular # fix Au angular momenta")

    return script

def generate_lammps_script(atoms, pctg, steady_state_steps, timestep):

    gold_atoms, gold_percentage, pd_percentage, daily_steps = compute_sim_features(atoms, pctg)
    t_initial = temp_dicts[atoms][pctg*100]

    dir_path = f"../inputs/set-8/mlmd/{round(timestep * 1000)}-fs/{atoms}-atoms/Pd{pd_percentage}Au{gold_percentage}"
    os.makedirs(dir_path, exist_ok=True)  # Create the directory if it doesn't exist

    daily_steps = steps_per_day[atoms]
    num_parts = -(-steady_state_steps // daily_steps)

    for part in range(1, num_parts + 1):
    
        previous_steps_taken = 0

        if part == 1:
            template_file = 'template_ss_mlmd.txt'
        else:
            template_file = 'template_restart_ss.txt'

        if part < num_parts:
            ss_time_part = daily_steps * (0.001 * timestep) # doing some ns to fs conversion. Hardcoded value is not the timestep, just converting multiplier 
        else:
            # Handle the last part specially if it doesn't align perfectly
            remaining_ss_steps = steady_state_steps - previous_steps_taken
            ss_time_part = remaining_ss_steps * (0.001 * timestep) # doing some ns to fs conversion. Hardcoded value is not the timestep, just converting multiplier 

        with open(template_file, 'r') as file:
            template = file.read()
            template = template.replace("${PD_PCTG}", str(pd_percentage))
            template = template.replace("${AU_PCTG}", str(gold_percentage))
            template = template.replace("${N_ATOMS}", str(atoms)) # Total number of atoms, different from gold_atoms which is only the number of gold atoms
            template = template.replace("${GOLD_ATOMS}", str(gold_atoms))

            template = template.replace("${T_T}", str(t_initial))

            template = template.replace("${TIMESTEP}", str(timestep))
            template = template.replace("${DUMP_FREQ}", str(round(1 / timestep)))
            template = template.replace("${DAMPING_FREQ}", str(round(timestep * 100, 1)))
            template = template.replace("${UNITARY_TIMESTEP}", str(round(timestep * 1000)))

            template = template.replace("${SS_TIME}", str(round(ss_time_part, 2)))
            template = template.replace("${SS_STEPS}", str(daily_steps if part < num_parts else remaining_ss_steps))
            template = template.replace("${RESTART_STEPS}", str(daily_steps if part < num_parts else remaining_ss_steps))
            

            template = template.replace("${N_SCRIPT}", str(part))
            template = momentum_fix(pctg, template)

            #print(template, "\n\n\n")

            output_name = f"prod_{atoms}_Pd{pd_percentage}Au{gold_percentage}_ss{'_restart' if part > 1 else ''}{'' if part == 1 else f'_{part - 1}'}.in"
            file_path = os.path.join(dir_path, output_name)

            with open(file_path, 'w') as output:
                output.write(template)

            print((f"File saved at: {file_path}"))
        
        previous_steps_taken += daily_steps if part < num_parts else remaining_ss_steps

for a in atoms:
    for p in percentages:
        generate_lammps_script(atoms = a, pctg = p, steady_state_steps = steady_state_steps, timestep = timestep)
import numpy as np
import os
from temp_dict import temps_55, temps_147, temps_309

atoms = [55, 147, 309]
percentages = [0.7, 0.76, 0.8, 0.86, 0.9, 0.96, 1]

timestep = 0.005 # in picosecond
t_floor = 300 

md_type = 'mlmd'

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

def nearest_multiple_20(n):
    if n % 20 == 0:
        return n
    else:
        difference = 20 - (n % 20)
        return n + difference
  
def compute_sim_features(atoms, pctg, t_floor):
    # Composition features
    gold_atoms = int(np.trunc(atoms * pctg))
    gold_percentage = int(round(pctg * 100))
    pd_percentage = 100 - gold_percentage
    gold_percentage = int(round(pctg * 100))

    temperature_dict = temp_dicts[atoms] 
    t_target = temperature_dict.get(gold_percentage)
    #t_target = int(round(nearest_multiple_20(t_target)))

    # Ramp features
    time_interval = (t_target - t_floor) / 20
    ramp_steps = int(round((time_interval) / (timestep * 0.001), 0))
    daily_steps = steps_per_day[atoms]

    return(gold_atoms, gold_percentage, pd_percentage, t_target, ramp_steps, daily_steps)

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

def generate_lammps_script(atoms, pctg, t_floor, timestep, md_type):
    gold_atoms, gold_percentage, pd_percentage, t_target, ramp_steps, daily_steps = compute_sim_features(atoms, pctg, t_floor)

    dir_path = f"../inputs/set-8/{md_type}/{round(timestep * 1000)}-fs/{atoms}-atoms/Pd{pd_percentage}Au{gold_percentage}"
    os.makedirs(dir_path, exist_ok=True)  # Create the directory if it doesn't exist

    # Calculate the number of parts the ramp needs to be split into
    num_parts = -(-ramp_steps // daily_steps)  # Ceiling division to handle any case of ramp_steps

    current_temperature = t_target
    previous_steps_taken = 0

    for part in range(1, num_parts + 1):
        template_file = 'template_restart_rampdown.txt'

        if part < num_parts:
            ramp_time_part = daily_steps * (0.001 * timestep) # doing some ns to fs conversion. Hardcoded value is not the timestep, just converting multiplier 
            temperature_decrease = ramp_time_part * 20
        else:
            # Handle the last part specially if it doesn't align perfectly
            remaining_ramp_steps = ramp_steps - previous_steps_taken
            ramp_time_part = remaining_ramp_steps * (0.001 * timestep) # doing some ns to fs conversion. Hardcoded value is not the timestep, just converting multiplier 
            temperature_decrease = ramp_time_part * 20

        t_final_part = int(np.round(current_temperature - temperature_decrease, 0))

        print(temperature_decrease)
        print(t_final_part)

        with open(template_file, 'r') as file:
            template = file.read()
            template = template.replace("${PD_PCTG}", str(pd_percentage))
            template = template.replace("${AU_PCTG}", str(gold_percentage))
            template = template.replace("${N_ATOMS}", str(atoms)) # Total number of atoms, different from gold_atoms which is only the number of gold atoms
            template = template.replace("${GOLD_ATOMS}", str(gold_atoms))

            template = template.replace("${T_FLOOR}", str(t_floor))
            template = template.replace("${T_IN}", str(current_temperature))
            template = template.replace("${T_F}", str(t_final_part))
            template = template.replace("${T_T}", str(t_target))

            template = template.replace("${TIMESTEP}", str(timestep))
            template = template.replace("${DUMP_FREQ}", str(round(1 / timestep)))
            template = template.replace("${DAMPING_FREQ}", str(round(timestep * 100, 1)))
            template = template.replace("${UNITARY_TIMESTEP}", str(round(timestep * 1000)))

            template = template.replace("${RAMP_DOWN_STEPS}", str(daily_steps if part < num_parts else remaining_ramp_steps))
            template = template.replace("${RAMP_DOWN_TIME}", str(round(ramp_time_part, 2)))
            template = template.replace("${RESTART_STEPS}", str(daily_steps if part < num_parts else remaining_ramp_steps))
            template = template.replace("${N_SCRIPT}", str(part))
            template = template.replace("${PICKUP}", "ss" if part == 1 else "rampdown")

            template = template.replace("${MD_TYPE}", str(md_type))

            template = momentum_fix(pctg, template)

            #print(template, "\n\n\n")

            output_name = f"prod_{atoms}_Pd{pd_percentage}Au{gold_percentage}_rampdown{'_restart' if part > 1 else ''}{'' if part == 1 else f'_{part - 1}'}.in"
            file_path = os.path.join(dir_path, output_name)

            with open(file_path, 'w') as output:
                output.write(template)

            print(f"File saved at: {file_path}")

        current_temperature = t_final_part
        previous_steps_taken += daily_steps if part < num_parts else remaining_ramp_steps

for a in atoms:
    for p in percentages:
        generate_lammps_script(atoms = a, pctg = p, t_floor = t_floor, timestep = timestep, md_type = md_type)
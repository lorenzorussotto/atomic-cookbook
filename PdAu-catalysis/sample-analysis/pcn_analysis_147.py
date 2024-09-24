import matplotlib
matplotlib.use('Agg')  # Set the backend to 'Agg' before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
import csv

# Load trajectory files, define simulation parameters
file_names = ['Pd30Au70_147_cuboct_730K_ss.all.bin', 
              "Pd24Au76_147_cuboct_720K_ss.all.bin",
              "Pd20Au80_147_cuboct_710K_ss.all.bin",
              "Pd14Au86_147_cuboct_700K_ss.all.bin",
              "Pd10Au90_147_cuboct_690K_ss.all.bin",
              "Pd04Au96_147_cuboct_670K_ss.all.bin"]

cutoff_distance = 3.5  # appropriate cutoff distance for metal-metal interactions in Pd-Au system

# Storage for time series data
results = {}

# Function to compute average coordination numbers
def compute_avg_coordination(data, cutoff, type_A, type_B):
    total_coordination = 0
    count_of_type_A = 0  # Track the number of type_A atoms

    finder = CutoffNeighborFinder(cutoff, data)
    for index in range(data.particles.count):
        if data.particles.particle_types[index] == type_A:
            count_of_type_A += 1
            total_coordination += sum(1 for n in finder.find(index) if data.particles.particle_types[n.index] == type_B)

    # Calculate average coordination if there are any atoms of type A
    if count_of_type_A > 0:
        avg_coordination = total_coordination / count_of_type_A
    else:
        avg_coordination = 0
    
    return avg_coordination

# Process trajectories
for file_name in file_names:

    pd_percentage = int(file_name[2:4])
    # Load simulation data
    print(f'Currently analyzing file {file_name}')
    pipeline = import_file(file_name)

    csv_data = []

    total_frames = pipeline.source.num_frames
    avg_Au_Au = avg_Pd_Pd = avg_Au_Pd = avg_Pd_Au = 0

    for frame in range(total_frames):
        
        data = pipeline.compute(frame)
        # Calculate average coordination numbers for each pair
        avg_Au_Au_frame = compute_avg_coordination(data, cutoff_distance, 2, 2)
        avg_Pd_Pd_frame = compute_avg_coordination(data, cutoff_distance, 1, 1)
        avg_Au_Pd_frame = compute_avg_coordination(data, cutoff_distance, 2, 1)
        avg_Pd_Au_frame = compute_avg_coordination(data, cutoff_distance, 1, 2)

        avg_Au_Au += avg_Au_Au_frame
        avg_Pd_Pd += avg_Pd_Pd_frame
        avg_Au_Pd += avg_Au_Pd_frame
        avg_Pd_Au += avg_Pd_Au_frame

        csv_data.append([avg_Au_Au_frame, avg_Pd_Pd_frame, avg_Au_Pd_frame, avg_Pd_Au_frame])

        if frame % 100 == 0:
            print(f'Analyzing Frame {frame}, Time: {round(frame / 1000, 4)} ns')
        
    # Average over all frames
    avg_Au_Au /= total_frames
    avg_Pd_Pd /= total_frames
    avg_Au_Pd /= total_frames
    avg_Pd_Au /= total_frames

    # Store results
    if pd_percentage not in results:
        results[pd_percentage] = []
    results[pd_percentage].append((avg_Au_Au, avg_Pd_Pd, avg_Au_Pd, avg_Pd_Au))

    # Saving for each composition
    csv_file_name = file_name.split('_ss')[0] + '.csv'  
    with open(csv_file_name, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['avg_Au_Au', 'avg_Pd_Pd', 'avg_Au_Pd', 'avg_Pd_Au']) # Header
        writer.writerows(csv_data)
        print(f'\n\nCSV file {csv_file_name} saved')

    print(f"\n\nTrajectory analysis for {file_name} is done")

pd_percentages = sorted(results.keys())
averages = {label: [] for label in ['Au-Au', 'Pd-Pd', 'Au-Pd', 'Pd-Au']}
for pd in pd_percentages:
    for avg in results[pd]:
        averages['Au-Au'].append(avg[0])
        averages['Pd-Pd'].append(avg[1])
        averages['Au-Pd'].append(avg[2])
        averages['Pd-Au'].append(avg[3])

# Saving data
average_pCN = "147_pCN_averages_ss.csv"
with open(average_pCN, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Pd_pctg'] + list(averages.keys()))
    for i, pd in enumerate(pd_percentages):
        row = [pd] + [averages[label][i] for label in averages]
        writer.writerow(row)
print(f"\n\nPlotted data saved to {average_pCN}")

# Plotting results
colors = {
    'Au-Au': 'darkorange',
    'Pd-Pd': 'darkblue',
    'Au-Pd': 'darkgreen',
    'Pd-Au': 'darkred'
}

plt.rcParams['savefig.dpi'] = 300  # High resolution
plt.figure(figsize=(10, 8))
for label, values in averages.items():
    plt.plot(pd_percentages, values, "-o", label=label, color=colors[label])
plt.xticks([4, 10, 14, 20, 24, 30])
plt.xlabel(r'$\mathrm{Pd}_x\mathrm{Au}_{1-x}$ (%)') 
plt.ylabel('Average p(CN)')
plt.title('147 atoms - Average pCN vs. Pd %')
plt.legend()
plt.savefig('pCN_vs_PdPercentage_147atoms.jpg')

print('\n\nFigure saved, job completed--')

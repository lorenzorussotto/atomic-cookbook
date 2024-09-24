from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier, ClusterAnalysisModifier, ExpressionSelectionModifier
from ovito.data import *
import numpy as np
import pandas as pd

# Define trajectories with their respective files
trajectories = [
    ['Pd04Au96_147_cuboct_670K_ss.all.bin', 'Pd04Au96_147_cuboct_670-300K_rampdown_mlmd_1.all.bin'],
    ['Pd10Au90_147_cuboct_690K_ss.all.bin', 'Pd10Au90_147_cuboct_690-300K_rampdown_mlmd_1.all.bin'], 
    ['Pd14Au86_147_cuboct_700K_ss.all.bin', 'Pd14Au86_147_cuboct_700-300K_rampdown_mlmd_1.all.bin'],
    ['Pd20Au80_147_cuboct_710K_ss.all.bin', 'Pd20Au80_147_cuboct_710-300K_rampdown_mlmd_1.all.bin'],
    ['Pd24Au76_147_cuboct_720K_ss.all.bin', 'Pd24Au76_147_cuboct_720-300K_rampdown_mlmd_1.all.bin'],
    ['Pd30Au70_147_cuboct_730K_ss.all.bin', 'Pd30Au70_147_cuboct_730-300K_rampdown_mlmd_1.all.bin']
]

results = []

# Process each pair of files in the trajectories
for files in trajectories:
    
    dfs = []

    for file in files:

        data_per_file = []
        # Load the file into the pipeline
        pipeline = import_file(file)

        # Add a modifier to calculate coordination numbers for all atoms
        pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=3.5))

        # Loop through frames
        for frame_index in range(pipeline.source.num_frames):

            # Apply cluster analysis for Pd atoms only in the entire NP
            pipeline.modifiers.append(ExpressionSelectionModifier(expression="ParticleType == 1"))
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=3.0, sort_by_size=True, only_selected=True))
            data = pipeline.compute(frame_index)

            # Get Pd cluster info for entire NP (only Pd atoms)
            total_clusters = np.bincount(data.particles['Cluster'], minlength=2)
            monomers_total = np.count_nonzero(total_clusters[1:] == 1)
            dimers_total = np.count_nonzero(total_clusters[1:] == 2)
            trimers_total = np.count_nonzero(total_clusters[1:] == 3)
            larger_clusters_total = np.count_nonzero(total_clusters[1:] > 3)

            # Remove previous modifiers for clear surface analysis
            pipeline.modifiers.pop()
            pipeline.modifiers.pop()

            # Count all surface particles
            pipeline.modifiers.append(ExpressionSelectionModifier(expression="Coordination < 10"))
            data = pipeline.compute(frame_index)
            surface_mask = data.particles['Selection']
            total_surface_atoms = np.sum(surface_mask)
            pipeline.modifiers.pop()

            # Count Pd surface particles only
            pipeline.modifiers.append(ExpressionSelectionModifier(expression="Coordination < 10 && ParticleType == 1"))
            data = pipeline.compute(frame_index)
            surface_pd_mask = data.particles['Selection']
            total_surface_Pd_atoms = np.sum(surface_pd_mask)

            # Finally, analysis of Pd surface clusters
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=3.0, sort_by_size=True, only_selected=True))
            data = pipeline.compute(frame_index)
            surface_clusters = np.bincount(data.particles['Cluster'], minlength=2)
            monomers_surface = np.count_nonzero(surface_clusters[1:] == 1)
            dimers_surface = np.count_nonzero(surface_clusters[1:] == 2)
            trimers_surface = np.count_nonzero(surface_clusters[1:] == 3)
            larger_clusters_surface = np.count_nonzero(surface_clusters[1:] > 3)

            # Save data for this frame
            frame_data = {
                'Frame': frame_index,
                'surface_monomers': monomers_surface,
                'surface_dimmers': dimers_surface,
                'surface_trimers': trimers_surface,
                'surface_four_plus_mers': larger_clusters_surface,
                'surface_pd_atoms': total_surface_Pd_atoms,
                'surface_total_atoms': total_surface_atoms,
                'total_monomers': monomers_total,
                'total_dimers': dimers_total,
                'total_trimers': trimers_total,
                'total_four_plus_mers': larger_clusters_total
            }

            data_per_file.append(frame_data)
            print(f'Analyzing Frame {frame_index} of file {file}')
            print(f'{frame_data}\n\n')  # Print current frame's data

            # Clear modifiers for the next frame to reset the selection and analysis
            pipeline.modifiers.clear()
            pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=3.5))

        # Collect data for this pair of files (trajectory)
        df = pd.DataFrame(data_per_file)
        dfs.append(df)
    
    final_trajectory_data = pd.concat(dfs, ignore_index=True)
    final_trajectory_data['Frame'] = range(1, len(final_trajectory_data) + 1)

    file_name = files[0].split('_')[0] + '_' + files[0].split('_')[1] + '_ensemble_data_set5'  # Extracts "PdXXAuYY_size"
    file_name = f"{file_name}.csv"
    final_trajectory_data.to_csv(file_name, index=False)
    print(f"Saved combined data for trajectory {files[0]} and {files[1]} as {file_name}")

print('\nAnalysis Completed')

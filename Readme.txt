
Input:
What you need to cluster & segment is a Spectrum (N_D * N_T2 matrix) for each participant
Also, you need the D_array and T2_array containing the values (coordinates) of the D & T2 dimensions. It is necessary for plotting the spectrum.
For our project, data of each participant is placed under a second-level directory gathered in a working directory (e.g. Participant/001_XXX). The codes are written according to such settings. All data are saved as matlab file.  You can change the code according to your preference.

(1) Spectrum_Pathology_Correlation.m
This code calculaes the (Spearmann) correlation between each meshpoint on the spectrum with histopathologic scores (activity & chronicity index). For each correlation to calculate, there should be N_patient data-points.
The output would be 2 correlation-coefficient map. Corresponding p-value maps are also given. All saved in a .mat file.

(2) Assign_Compartments.m
This code assign the meshpoints into "compartment" (subregions on the spectrum). 
Please see our paper to find details.
The output would be a D_T2_cluster matrix (N_D * N_T2 * N_compartment, value: 0 or 1) and a D_T2_map matrix (N_D * N_T2, value: 1, 2, 3 ... N_compartment) 

(3) Batch_Segment.m
This code segements the spectrum into compartments, cyclying through all the participants
Only the spectrum is necessary input (even participant-level spectrum if it is the case) 
The output is an excel file containing volume fractions for each participant

(4) Batch_plot_Fmaps.m
This code plots the volume fraction maps of each compartment based on the segmentation strategy (D_T2_cluster matrix)
It cyclyes through all the participants
Necessary input: coordinate list of the voxels in ROI; voxel-level spectra within the ROI; segmentation strategy 
Figures would be saved under figure directory of each participant directory

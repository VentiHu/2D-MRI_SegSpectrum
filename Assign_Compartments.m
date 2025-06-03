close all
clear


%% Spectrum-Pathology Correlation-Based Segmentation
% This script defines the D_T2_LabelMap matrix based on spectrum-pathology correlations.

% Load the required data
load('Data_Spearman.mat'); 
% In this datafile, D_array, T2_array, rhoMatrix_AI, rhoMatrix_CI are loaded
% D_array should be an array containing the meshpoint values of diffusivity (D-axis on the spectrum)
% T2_array should be an array containing the meshpoint values of T2 (T2-axis on the spectrum)
% rhoMatrices are defined in spectrum_pathology_correlation.m

% Convert D_array to μm^2/ms as needed
D_array = D_array * 1000;

% Initialize the cluster matrix
N_D = length(D_array);
N_T2 = length(T2_array);
D_T2_LabelMap = zeros(N_D, N_T2);

% Define compartment identifiers
CompartmentNames = ["S", "R", "L", "PD"];
PD = 4; % Pseudo-diffusion compartment
S = 1;  % Short relaxation compartment
R = 2;  % Restricted diffusivity
L = 3;  % Long relaxation

% Loop through each mesh point and assign compartments
for i = 1:N_D
    for j = 1:N_T2
        D_value = D_array(i);
        T2_value = T2_array(j);
        rho_AI = rhoMatrix_AI(i, j);
        rho_CI = rhoMatrix_CI(i, j);
        
        % Assign Pseudo-diffusion compartment (PD)
        if (D_value > 4.0 && rho_CI < 0 && abs(rho_AI) <= abs(rho_CI)) || D_value > 15
            D_T2_LabelMap(i, j) = PD;
        % Assign Short relaxation compartment (S)
        elseif D_value < 10 && T2_value < 60 && rho_AI < 0 && rho_AI < rho_CI
            D_T2_LabelMap(i, j) = S;
        % Assign Long relaxation with free diffusion (L)
        elseif D_value < 10 && T2_value > 30 && rho_CI>0 && abs(rho_CI) > abs(rho_AI)
            D_T2_LabelMap(i, j) = L;
        % Assign Extended relaxation with restricted diffusivity (R)
        elseif D_value < 1.2 && T2_value > 20 && T2_value < 100 && rho_CI > 0 && abs(rho_AI) > abs(rho_CI)
            D_T2_LabelMap(i, j) = R;

        end
    end
end

% Region Growing for Unassigned Mesh Points
D_T2_LabelMap = regionGrowingForUnassigned(D_T2_LabelMap, N_D, N_T2);


% Transfer into D_T2_cluster
N_peak = max(max(D_T2_LabelMap));
D_T2_cluster = zeros(N_D,N_T2,N_peak);
for j = 1:N_peak  
    D_T2_cluster(:,:,j) = (D_T2_LabelMap == j);
end

% Visualize the D_T2_LabelMap matrix
figure;
imagesc(T2_array, D_array, D_T2_LabelMap);
xlabel('T2 (ms)', 'FontSize', 14);
ylabel('D (\mum^2/ms)', 'FontSize', 14);
set(gca,'FontSize',14,'XScale','log','YScale','log','YDir','Normal');
axis([min(T2_array) max(T2_array) min(D_array) max(D_array)]);
axis square
%colorbar;
colormap(jet);
caxis([0 4.5])
title('Spectral Region Assignment', 'FontSize', 16);

% Save the result
save('D_T2_cluster_result.mat', 'D_T2_LabelMap','D_T2_cluster', 'CompartmentNames');

%% Subfunction for Region Growing
function D_T2_LabelMap = regionGrowingForUnassigned(D_T2_LabelMap, N_D, N_T2)
    % Iterate through unassigned points and assign them based on neighboring regions
    unassigned_points = find(D_T2_LabelMap == 0);
    while ~isempty(unassigned_points)
        for idx = 1:length(unassigned_points)
            [i, j] = ind2sub(size(D_T2_LabelMap), unassigned_points(idx));
            neighbors = [];
            
            % Collect neighbors' indices (8-connectivity)
            if i > 1
                neighbors = [neighbors; i-1, j];
            end
            if i < N_D
                neighbors = [neighbors; i+1, j];
            end
            if j > 1
                neighbors = [neighbors; i, j-1];
            end
            if j < N_T2
                neighbors = [neighbors; i, j+1];
            end
            if i > 1 && j > 1
                neighbors = [neighbors; i-1, j-1];
            end
            if i > 1 && j < N_T2
                neighbors = [neighbors; i-1, j+1];
            end
            if i < N_D && j > 1
                neighbors = [neighbors; i+1, j-1];
            end
            if i < N_D && j < N_T2
                neighbors = [neighbors; i+1, j+1];
            end
            
            % Check the compartment labels of the neighbors
            neighbor_labels = [];
            for n = 1:size(neighbors, 1)
                ni = neighbors(n, 1);
                nj = neighbors(n, 2);
                if D_T2_LabelMap(ni, nj) ~= 0
                    neighbor_labels = [neighbor_labels; D_T2_LabelMap(ni, nj)];
                end
            end
            
            % Assign the current point to the most frequent neighbor label
            if ~isempty(neighbor_labels)
                D_T2_LabelMap(i, j) = mode(neighbor_labels);
            end
        end
        % Update the list of unassigned points
        unassigned_points = find(D_T2_LabelMap == 0);
    end
end


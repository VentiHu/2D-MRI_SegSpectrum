clear
close all

load('../D_T2_cluster_result.mat')
load('../Data_Spearman.mat');
% CompartmentNames could be loaded from D_T2_cluster.mat, or defined here
%CompartmentNames = ["S", "R", "L", "PD"];

% Convert D_array to μm^2/ms as needed
D_array = D_array * 1000;
N_D = length(D_array);
N_T2 = length(T2_array);
N_comp = size(D_T2_cluster,3);

% Create D and T2 matrix, representing D and T2 value of each mesh-point
D_matrix = repmat(D_array,1,N_T2);
T2_matrix = repmat(T2_array',N_D,1);

for c = 1:N_comp
    N_meshpoints = sum(sum(D_T2_cluster(:,:,c)));
    disp("Compartment "+CompartmentNames(c)+" N_meshpoints: "+num2str(N_meshpoints))
    
    D_collect = D_matrix .* D_T2_cluster(:,:,c);
    D_list = nonzeros(D_collect);
    %D_geo_ave = prod(D_list)^(1/length(D_list));
    D_geomean = geomean(D_list);
    disp("Compartment "+CompartmentNames(c)+" D_geomean: "+num2str(D_geomean))
    
    T2_collect = T2_matrix .* D_T2_cluster(:,:,c);
    T2_list = nonzeros(T2_collect);
    %T2_geo_ave = prod(T2_list)^(1/length(T2_list));
    T2_geomean = geomean(T2_list);
    disp("Compartment "+CompartmentNames(c)+" T2_geomean: "+num2str(T2_geomean))
    
end

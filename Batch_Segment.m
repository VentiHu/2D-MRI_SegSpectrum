clear
close all
clc

outputDir = pwd();
workdir = "..\..\Participants"; % the directory saving te patient files

load('D_T2_cluster_result.mat') % Created by "Assign_Compartment.m", or from your own segmentation strategy
CompartmentNames = ["S", "R", "L", "PD"]; % This could also be loaded from D_T2_cluster_result.mat

cd(workdir)
% Catch all the participant datafile directories
subdir = dir(pwd());
num_patient = 0;
name_patient = strings(0);
for i = 1:length(subdir)
	if(	isequal(subdir(i).name, '.' )||...
		isequal(subdir(i).name, '..')||...
		~subdir(i).isdir)
		continue;
	end
	num_patient = num_patient + 1;
	name_patient = [name_patient,subdir(i).name];
end

num_patient
mainfig_dir = "figure";
N_peak = 4; % Number of compartments
list_ROI_ratio = zeros(num_patient,N_peak); %patient-list of all compartments ratios
list_ROI_SliceRatio = zeros(num_patient,N_peak);

i_slice = 0;
name_slice = strings(num_patient,1);

%OA_info = readmatrix('OA.xlsx');

%% MAIN LOOP
for n = 1:num_patient
	disp(name_patient(n))
	cd(name_patient(n))
    clear Spectrum_voxel_norm
    
    load('ROI.mat')
    % Contains: mask_ROI (3D matrix)
    load('Step1_data.mat') 
    % Contains: img_b0 (DWI b0 image)
    % Contains: ROI_points = number of voxels in ROI
    % Contains: ROI_coor = coordinate of voxels in ROI; ROI_points * 3 matrix
    load('Step2_data.mat') 
    % Contains: Spectrum_voxel_norm_sum (N_D * N_T2 * N_ROIvoxel)
    % Contains: D_array (N_D * 1 array, output of MERA tool)
    % Contains: T2_array (N_T2 * 1 array, output pf MERA tool)
    
    [len_x, len_y, ~] = size(img_b0);
    Spectrum_voxel_norm_sum = Spectrum_voxel_norm_sum/sum(sum(Spectrum_voxel_norm_sum)); % Normalization
    
    %% SEGEMENTATION PARAMETER
    N_D = size(D_array,1);
    N_T2 = size(T2_array,1);	

    % OBTAIN OrderMatrix OF WHOLE KIDNEY ROI
	% FOR EXAMPLE,
	% 0 0 0 1 2 0 0
	% 0 0 3 4 5 6 0
	% 0 0 7 8 9 0 0
	% ......
	OrderMatrix_ROI = zeros(len_x, len_y, num_slice);
	for j = 1:ROI_points
		x = ROI_coor(j,1);
		y = ROI_coor(j,2);
		z = ROI_coor(j,3);
		OrderMatrix_ROI(x,y,z) = j;
    end

	%% SEGMENT THE SUB-REGIONAL SPECTRA
    list_ROI_ratio(n,:) = mySegment(Spectrum_voxel_norm_sum,D_T2_cluster);
 
    %% PLOT SPECTRA
    my_plot(Spectrum_voxel_norm_sum,D_array,T2_array,mainfig_dir+"\Spectrum.png")
    
    %% SLICE LOOP
    %{
    for zi = 1:num_slice
        if sum(sum(mask_ROI(:,:,zi)))==0
            continue;
        end
        i_slice = i_slice + 1;
        name_slice(i_slice) = name_patient(n) + " slice"+ num2str(zi);
        slice_prefix = mainfig_dir + "\slice" + num2str(zi);
        
        
        % CALCULATE SLICE AVERAGE
        mask_ROI_ThisSlice = zeros(len_x,len_y,num_slice);
        mask_ROI_ThisSlice(:,:,zi) = mask_ROI(:,:,zi);
        Spectrum_ROI_ThisSlice = mySubregionSpectrum(mask_ROI_ThisSlice,OrderMatrix_ROI,Spectrum_voxel_norm);
        list_ROI_SliceRatio(i_slice,:) = mySegment(Spectrum_ROI_ThisSlice,D_T2_cluster);
        
    end
    %}
    close all
    cd ..
end

cd(outputDir)

varNames = ["Name","ROI_VF","cortex_VF","medulla_VF",];
table_results = table(name_patient.',list_ROI_ratio,list_cortex_ratio,list_medulla_ratio,'VariableNames',varNames);
writetable(table_results,'results.xlsx','sheet','patients');

% Slice-level output
%{
table_results_slice = table(name_slice,list_ROI_SliceRatio,list_cortex_SliceRatio,list_medulla_SliceRatio,'VariableNames',varNames);
writetable(table_results_slice,'results.xlsx','sheet','slices');
%}
disp('Finish')



%% FUNCTION DEFINITION -----------------------------------------------
function [fraction] = mySegment(spectrum,matrix_div_spec)
    % THIS FUNCTION CALCULATE THE VOLUME FRACTION OF EACH COMPARTMENT BY 
    % DEFINED SPECTRUM DIVISION MATRIX
	[N_mesh1,N_mesh2,N_peak] = size(matrix_div_spec);
	sumI_peak = zeros(1,N_peak); % I = INTENSITY
	fraction = zeros(1,N_peak);
	
	for j = 1:N_peak
		sumI_peak(j) = sum(sum(spectrum .* matrix_div_spec(:,:,j)));
    end
    fraction = sumI_peak * 100;
end

function [] = my_plot(PlotSpectrum,D_array,T2_array,save_name)
	figure
	Maxvalue = max(max(PlotSpectrum));
    Minvalue = min(min(PlotSpectrum));
    Boundary = max(Maxvalue,abs(Minvalue));
    hcon = contour(T2_array,D_array*1000,PlotSpectrum,20,'linewidth',2);
	xlabel('T2 (ms)','FontSize',14)
	ylabel('D (10^-^3 mm^2/s)','FontSize',14)
	grid on, view([0 90]), shading interp
	set(gca,'FontSize',14,'XScale','log','YScale','log');
	axis([min(T2_array) max(T2_array) min(D_array*1000) max(D_array*1000)]);
	colorbar;
	saveas(gcf,save_name)
	close gcf
end

function [Spectrum_output] = mySubregionSpectrum(mask,OrderMatrix_ROIROI,Spectrum_ROI)
    % THIS FUNCTION CALCULATE THE SUB-REGIONAL VOXEL-SUM SPECTRUM
    [N_mesh1,N_mesh2,~] = size(Spectrum_ROI);
	OrderMatrix = mask .* OrderMatrix_ROIROI;
	OrderArray = OrderMatrix(:);
	OrderArray(OrderArray==0)=[];
	OrderArray = sort(OrderArray)';
	Spectrum_output = zeros(N_mesh1,N_mesh2);
	for j = OrderArray
		Spectrum_output = Spectrum_output + Spectrum_ROI(:,:,j);
    end
    Spectrum_output = Spectrum_output / length(OrderArray);
end
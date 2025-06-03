clear
close all
clc

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
% 
fig_dir = "figure"
N_peak = 4;

%% MAIN LOOP
for n = 1:num_patient
	disp(name_patient(n))
	cd(name_patient(n))
%	rmdir('_figure2','s')
	mkdir(fig_dir)
	
	%% data preparation
	load('ROI.mat');
    % Contains: mask_ROI (3D matrix)
	load('Step1_data.mat'); 
    % Contains: img_b0 (DWI b0 image)
    % Contains: ROI_points = number of voxels in ROI
    % Contains: ROI_coor = coordinate of voxels in ROI; ROI_points * 3 matrix
	load('Step2_data.mat');
	% Contains: Spectrum_voxel_norm_sum (N_D * N_T2 * N_ROIvoxel matrix)
    
	[len_x, len_y, ~] = size(img_b0);
	b_array = [0 50 100 300 600 1200];
	TE_array = [55 70 95 120 145 170];

	N_D = size(D_array,1);
	N_T2 = size(T2_array,1);
	img_black_plot = zeros(len_x,len_y,3);
	[xx, yy] = meshgrid(1:len_x, 1:len_y);

    for zi = 1:num_slice % Slice Loop
        if sum(sum(mask_ROI(:,:,zi))) == 0
            continue; % PLOT IMAGES ONLY FOR SLICES WITH ROI AREA > 0
        end
        slice_dir = "slice" + num2str(zi);
        slice_prefix = fig_dir + "\slice" + num2str(zi);
        if ~exist(slice_dir,'dir')==0
            rmdir(slice_dir,'s')
        end
        %mkdir(slice_dir)
        %% PLOT DWI & ROI

        figure
        imshow(img_b0(:,:,zi),[]);
        saveas(gcf,slice_prefix+"_DWI_b0.png")

        figure
        imshow(mask_ROI(:,:,zi),[]);
        saveas(gcf,slice_prefix+"_mask.png")

    %}	


        %% Plot Compartments

        spec_voxel_peak = zeros(N_D,N_T2,ROI_points,N_peak); %single-voxel spectrum segmented into each peak
        sum_spec_voxel_peak = zeros(ROI_points,N_peak);
        frac_voxel_peak = zeros(ROI_points,N_peak); %fractions of each peak in single voxel
        Vmap_plot = zeros(len_x,len_y,N_peak); %map matrix

        for i = 1:ROI_points % ROI_points = number of voxels in ROI
            if ROI_coor(i,3) == zi %  ROI_coor = coordinate of voxels in ROI; ROI_points * 3 matrix
                % COMPARTMENT SUM UP
                for j = 1:N_peak
                    spec_voxel_peak(:,:,i,j) = Spectrum_voxel(:,:,i).*D_T2_cluster(:,:,j);
                    sum_spec_voxel_peak(i,j) = sum(sum(spec_voxel_peak(:,:,i,j)));
                end
                % CALCULATE MAPPING
                if sum(sum_spec_voxel_peak(i,:),2) == 0
                    frac_voxel_peak(i,j) = 0;
                else
                    for j = 1:N_peak
                        frac_voxel_peak(i,j) = sum_spec_voxel_peak(i,j)/sum(sum(Spectrum_voxel(:,:,i)));
                        if frac_voxel_peak(i,j) == nan
                            frac_voxel_peak(i,j) = 0;
                        end
                        xi = ROI_coor(i,1);
                        yi = ROI_coor(i,2);
                        Vmap_plot(xi,yi,j) = frac_voxel_peak(i,j);
                    end
                end
            end
        end
        % Create background from DWI b0 figures
        % Alternatively, you can use img_black_plot as background
        img_b0_plot = zeros(len_x,len_y,3);
        img_b0_double = im2double(img_b0(:,:,zi));
        maxValue_ROI = max(max(img_b0_double));
        for j = 1:3
            img_b0_plot(:,:,j) = img_b0_double/maxValue_ROI;
        end
        
        % PLOT
        map_to_plot = zeros(len_x,len_y);
        for j = 1:N_peak
            map_to_plot = Vmap_plot(:,:,j);
            map_to_plot(find((map_to_plot+mask_ROI(:,:,zi))==0)) = nan;
            
            my_plot_VFmap(map_to_plot,"VF_"+CompartmentNames(j),img_b0_plot)
            
            str_FigSave = slice_prefix+"_VF_"+CompartmentNames(j)+".png";
            saveas(gcf,str_FigSave)
        end
    end

	cd ..
	close all
end

disp('Finish')

% SUBFUNCTION
function output = my_plot_VFmap(VFmap,fig_name,img_background)
    
    [len_x, len_y] = size(VFmap);
    [xx, yy] = meshgrid(1:len_x, 1:len_y);
    
    figure('Name',fig_name)
    imshow(img_background,[]);
    hold on;
    pcolor(yy,xx,VFmap'*100);
    shading interp
    colormap(jet)
    colorbar;
    set(gca,'position',[0.1 0.1 0.8 0.8]);
    title(fig_name);
    set(gca,'YDir','reverse');
    caxis([0 100]);

end
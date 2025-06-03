clear
close all
clc

outputDir = pwd();
workDir = "..\..\Participants"; % the directory saving te patient files
cd(workDir)

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
num_inclusion = 66; % include only LN, no HV

%% PATHOLOGY PREPARATION
info = readtable('patho.xlsx','sheet',1); % Table containing the pathology info, columns: AI, CI
% AI = Activity Index; CI = Chronicity Index
ActiveIndex = patho.AI;
ChronicIndex = patho.CI;


%% SPECTRUM PREPARATION
N_D = 30; % D = Diffusivity
N_T2 = 30; % T2 = T2 relaxation time
Spectrum_All = zeros(N_D,N_T2,num_inclusion);
Spectrum_temp = zeros(N_D,N_T2);

for n = 1:num_inclusion
	disp(name_patient(n))
	cd(name_patient(n))
		

    load('Step2_data.mat')
    % Contains: Spectrum_voxel_norm_sum (N_D * N_T2 * N_ROIvoxel matrix)
    Spectrum_temp = Spectrum_voxel_norm_sum;
    
    % RENORMALIZE SUM TO 1
	Amplitude = sum(sum(Spectrum_temp));
    Spectrum_All(:,:,n) = Spectrum_temp / Amplitude;

	cd ..
end

Spectrum_Ave = mean(Spectrum_All,3);
Spectrum_Ave = Spectrum_Ave / max(max(Spectrum_Ave));

%% Do Spectrum-Pathology Spearman Correlation
% AI = Activity Index; CI = Chronicity Index
[rhoMatrix_AI, pvalMatrix_AI] = do_spectral_Spearman(Spectrum_All, ActiveIndex);
[rhoMatrix_CI, pvalMatrix_CI] = do_spectral_Spearman(Spectrum_All, ChronicIndex);

%% PLOT
figure;
h = heatmap(flipud(Spectrum_Ave));
h.Title = 'Subject-Averaged Spectrum';
h.XLabel = 'T2';
h.YLabel = 'Diffusivity';
h.Colormap = hot;
h.ColorLimits = [0, 1]; 
saveas(h,outputDir+"Spectrum_All_heatmap-style")

% Plot correlation heatmap & p-value mapping

% Active Index
my_plot_heatmap(rhoMatrix_AI,D_array,T2_array,outputDir+"\Rho_AI_heatmap.png")
my_plot_pvalmap(pvalMatrix_AI,D_array,T2_array,outputDir+"\P-val_AI_heatmap.png")
% Chronic Index
my_plot_heatmap(rhoMatrix_CI,D_array,T2_array,outputDir+"\Rho_CI_heatmap.png")
my_plot_pvalmap(pvalMatrix_CI,D_array,T2_array,outputDir+"\P-val_CI_heatmap.png")

%
% Save the statistic results
save('Data_Spearman.mat','Spectrum_All','Spectrum_Ave','rhoMatrix_AI','pvalMatrix_AI','rhoMatrix_CI','pvalMatrix_CI','D_array','T2_array')

%% SUBFUNCTIONS --------------------------------------
function [rhoMatrix, pvalMatrix] = do_spectral_Spearman(spectrum_series, pathology)
% Input:
% Suppose we have N subjects, each with a sz1 * sz2 spectrum
% spectrum_series is a sz1 * sz2 * N matrix
% pathoglogy is an array with length of N

[sz1, sz2, N] = size(spectrum_series);
rhoMatrix = zeros(sz1,sz2);
pvalMatrix = zeros(sz1,sz2);
    for i = 1:sz1
        for j = 1:sz2
            data = reshape(spectrum_series(i,j,:),[N,1]);
            % Calculate the Spearman correlation coefficient
            [rho, pval] = corr(data, pathology, 'Type', 'Spearman');
            rhoMatrix(i,j) = rho;
            pvalMatrix(i,j) = pval;
        end
    end
end

function output = my_plot_heatmap(PlotMatrix,D_array,T2_array,save_name)
    % DEFINE COLORMAP
    mycolorpoint=[[0 0 255];[255 255 255];[255 0 0]];
    mycolorposition=[1 32 63];
    mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:63,'linear','extrap');
    mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:63,'linear','extrap');
    mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:63,'linear','extrap');
    mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/255;
    color_heatmap = round(mycolor*10^4)/10^4;

    figure; 
    h = imagesc(T2_array,D_array*1000,PlotMatrix);
	xlabel('T2 (ms)','FontSize',14)
	ylabel('D (10^-^3 mm^2/s)','FontSize',14)
	set(gca,'FontSize',14,'XScale','log','YScale','log','YDir','Normal');
	axis([min(T2_array) max(T2_array) min(D_array*1000) max(D_array*1000)]);
    colorbar
    colormap(color_heatmap);
    caxis([-1 1]); % Sets the color limits for the heatmap
    saveas(h,save_name)
end

function output = my_plot_pvalmap(PlotMatrix,D_array,T2_array,save_name)
    % DEFINE COLORMAP    
    mycolorpoint=[[0 255 0];[220 255 220];[255 255 255]];
    mycolorposition=[1 5 63];
    mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:63,'linear','extrap');
    mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:63,'linear','extrap');
    mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:63,'linear','extrap');
    mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/255;
    color_pval = round(mycolor*10^4)/10^4;
    
    figure; 
    h = imagesc(T2_array,D_array*1000,PlotMatrix);
	xlabel('T2 (ms)','FontSize',14)
	ylabel('D (10^-^3 mm^2/s)','FontSize',14)
    set(gca,'FontSize',14,'XScale','log','YScale','log','YDir','Normal');
	axis([min(T2_array) max(T2_array) min(D_array*1000) max(D_array*1000)]);
    colorbar
    colormap(color_pval);
    caxis([0 1]); % Sets the color limits for the heatmap
    saveas(h,save_name)
end

clear; close all; clc

% Parameters
nBootstraps = 100;  % Number of bootstrap samples
outputDir = pwd();
workDir = "..\..\..\BatchProcess_Modified";
cd(workDir)

% Load patient info
subdir = dir(pwd()); num_patient = 0; name_patient = strings(0);
for i = 1:length(subdir)
    if ~subdir(i).isdir || isequal(subdir(i).name, '.') || isequal(subdir(i).name, '..')
        continue;
    end
    num_patient = num_patient + 1;
    name_patient(end+1) = subdir(i).name;
end

% Patient exclusion settings
num_inclusion = 66; 
list_exclude = [38];  % Excluded patients
list_nnls = [32,51,60,66];

%% PATHOLOGY PREPARATION
info = readtable('info.xlsx','sheet',1);
ActiveIndex = info.AI(1:num_inclusion); ActiveIndex(list_exclude) = [];
ChronicIndex = info.CI(1:num_inclusion); ChronicIndex(list_exclude) = [];

%% SPECTRUM PREPARATION
N_ADC = 30; N_T2 = 30;
Spectrum_All = zeros(N_ADC, N_T2, num_inclusion);
Spectrum_temp = zeros(N_ADC, N_T2);

% Load spectral data
for n = 1:num_inclusion
    cd(name_patient(n))
    if ismember(n, list_nnls)
        load('Step2_data_nnls.mat');
        Spectrum_temp = Spectrum_voxel_norm_sum;
    else
        load('Step2_data.mat');
        Spectrum_temp = Spectrum_voxel_norm_sum;
    end
    Amplitude = sum(sum(Spectrum_temp));
    Spectrum_All(:,:,n) = Spectrum_temp / Amplitude;
    cd ..
end
Spectrum_All(:,:,list_exclude) = [];  % Remove excluded patients

% Prepare grid coordinates
%[D_array, T2_array] = prepareGridCoordinates();  % Implement your grid definition

%% BOOTSTRAP RESAMPLING
% Initialize storage for bootstrap results
rhoMatrix_AI_all = zeros(N_ADC, N_T2, nBootstraps);
rhoMatrix_CI_all = zeros(N_ADC, N_T2, nBootstraps);

% Main bootstrap loop
for boot_iter = 1:nBootstraps
    % 1. Resample patients with replacement
    boot_idx = datasample(1:size(Spectrum_All,3), size(Spectrum_All,3), 'Replace', true);
    
    % 2. Extract resampled data
    boot_spectrum = Spectrum_All(:,:,boot_idx);
    boot_AI = ActiveIndex(boot_idx);
    boot_CI = ChronicIndex(boot_idx);
    
    % 3. Compute correlation matrices for this bootstrap sample
    [rhoMatrix_AI_all(:,:,boot_iter), ~] = do_spectral_Spearman(boot_spectrum, boot_AI);
    [rhoMatrix_CI_all(:,:,boot_iter), ~] = do_spectral_Spearman(boot_spectrum, boot_CI);
    
    fprintf('Bootstrap iteration %d/%d completed\n', boot_iter, nBootstraps);
end

%% AVERAGE RESULTS
% Calculate mean correlation matrices
rhoMatrix_AI_avg = mean(rhoMatrix_AI_all, 3);
rhoMatrix_CI_avg = mean(rhoMatrix_CI_all, 3);

%% PLOT AVERAGED RESULTS

% Plot averaged correlation maps
my_plot_heatmap(rhoMatrix_AI_avg, D_array, T2_array, fullfile(outputDir, 'Rho_AI_bootstrap_avg.png'));
my_plot_heatmap(rhoMatrix_CI_avg, D_array, T2_array, fullfile(outputDir, 'Rho_CI_bootstrap_avg.png'));

% Save bootstrap results
save(fullfile(outputDir, 'Data_Spearman_Bootstrap.mat'), ...
    'rhoMatrix_AI_all', 'rhoMatrix_CI_all', ...
    'rhoMatrix_AI_avg', 'rhoMatrix_CI_avg', ...
    'D_array', 'T2_array');

%% SUBFUNCTIONS
% ... (Keep existing subfunctions unchanged)

function [rhoMatrix, pvalMatrix] = do_spectral_Spearman(spectrum_series, pathology)
    [sz1, sz2, N] = size(spectrum_series);
    rhoMatrix = zeros(sz1, sz2);
    pvalMatrix = zeros(sz1, sz2);
    
    for i = 1:sz1
        for j = 1:sz2
            data = reshape(spectrum_series(i,j,:), [N,1]);
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
    color_heatmap = round(mycolor*10^4)/10^4;%保留4位小数

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
    color_pval = round(mycolor*10^4)/10^4;%保留4位小数
    
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
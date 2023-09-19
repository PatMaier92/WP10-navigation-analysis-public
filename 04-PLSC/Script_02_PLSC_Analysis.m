clear; close all; clc; format long;
addpath(genpath(pwd)) % add subfolder functions to path

%% Partial Least Square Correlation (PLSC) Analysis 
% @ date 2022-09-05 @ author Patrizia Maier & tracked by git 

% This script applies Partial Least Squares Correlation (PLSC; Krishnan et al., 2011) 
% to extract a latent variable (LV) capturing age-related differences in
% navigation behavior and correlating this LV with memory. 

% This script requires the plscmd toolbox (https://www.rotman-baycrest.on.ca).

%%  Load data and run PLSC 
%--------------------------------------------------------------------------
% LOAD DATA 
%--------------------------------------------------------------------------
% load raw data 
path = '../WP10_data/WP10_results/';

% age by session 1 learning trials
data_age_by_NlS1 = readtable([path, 'wp10_plsc_by_age.txt']); 

% memory data 
memory_table = readtable([path, 'wp10_plsc_memory.txt']); 

% settings 
data_cell = { data_age_by_NlS1 };
data_names = { 'by_age' }; 
analysis_method = { 3 }; % 3=regular behavior PLS; 5=non-rotated behavior PLS
n_conditions = { 1 }; % number of conditions 
by_group = { 0 }; % 0=no; 1=yes; note: not required if age is outcome variable 
design_contrasts = { [] }; % rows=conditions x groups, sorted as condition in group; colums = desired contrasts

for i=1:numel(data_cell)
    %--------------------------------------------------------------------------
    % CONFIGURATION
    %-------------------------------------------------------------------------- 
    data = data_cell{i};
    file_name = data_names{i};
    
    % behavioral output data
    plsinput.y = data.age;
    
    % behavioral explanatory data
    plsinput.X = table2array(data(:,5:size(data,2)));

    % z-standardization
    plsinput.y = zscore(plsinput.y,0,1);
    plsinput.X = zscore(plsinput.X,0,1);

    % cfg settings
    cfg.pls = [];
    cfg.pls.num_perm = 5000; % number of permutations
    cfg.pls.num_boot = 5000; % number of bootstrap tests
    cfg.pls.clim     = 95; % confidence interval level

    % analysis method: 3=regular behavior PLS; 5=non-rotated behavior PLS
    cfg.pls.method   = analysis_method{i}; 
    
    % add behavioral output data
    cfg.pls.stacked_behavdata = plsinput.y;
    
    % number of conditions 
    n_con = n_conditions{i}; 
    if cfg.pls.method==5
        cfg.pls.stacked_designdata=design_contrasts{i};
    end 

    % number of subjects and data preparation
    if by_group{i}==0
        n_subj = size(plsinput.y,1) / n_con;
        datamat1_allgroups = plsinput.X;
    elseif by_group{i}==1
        n_subj = histc(data(:,2),unique(data(:,2))) / n_con;
        datamat1_group1 = plsinput.X(1:n_subj(1),:);
        datamat1_group2 = plsinput.X(n_subj(1)+1:n_subj(1)+n_subj(2),:);
        datamat1_group3 = plsinput.X(n_subj(1)+n_subj(2)+1:end,:);
    end   
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % RUN PLS AND SAVE OUTPUT
    %--------------------------------------------------------------------------
    % input arguments: data, number of subjects, number of conditions, specific settings
    if by_group{i}==0
        plsres = pls_analysis({ datamat1_allgroups }, n_subj, n_con, cfg.pls);
    elseif by_group{i}==1
        plsres = pls_analysis({ datamat1_group1, datamat1_group2, datamat1_group3 }, n_subj, n_con, cfg.pls);
    end
    
    % add relevant input data to pls output 
    plsres.data.group = data.group; 
    plsres.data.age = data.age; 
    plsres.data.memoryAvg = memory_table.memoryAvg; 
    plsres.data.memoryEgo1 = memory_table.memoryEgo1; 
    plsres.data.memoryEgo2 = memory_table.memoryEgo2; 
    plsres.data.memoryAllo1 = memory_table.memoryAllo1; 
    plsres.data.memoryAllo2 = memory_table.memoryAllo2; 

    % calculate latent profile score (usc)
    plsres.usc_nav = datamat1_allgroups(:,1:4) * plsres.u(1:4); 
       
    % save data file 
    % save([path, '/wp10_plsc_', file_name, '_results_m', int2str(cfg.pls.method), '_g', int2str(size(n_subj,1)), '.mat'],'plsres');
    save([path, '/wp10_plsc_', file_name, '.mat'], 'plsres');
    
    clear n_con cfg datamat* p n_subj data file_name plsinput plsres;
end

clearvars; 
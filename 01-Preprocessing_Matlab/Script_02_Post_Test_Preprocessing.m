clear; close all; clc; format long;
addpath(genpath(pwd)) % add subfolder functions to path

%% PostTest Data Processing
% @ date 2021-05-21 @ author Patrizia Maier & tracked by git 
% for starmaze version WP10 Frankfurt

% The script requires .csv files as input.
% One trial_results file with general information (id, session, condition).

% Block 1: Set up input/output folders, Starmaze and Practise environment
% Block 2: Data analysis
% Block 3: Write data to xlsx file

%% Block 1: Set up input/output folders
%% data input folder and participant information
[data_folder] = setInputPath(); % provide folder with all raw data
[participant_start,participant_end] = setParticipants(); % provide participant range
s=4;

%% result folder
result_folder=[data_folder '\WP10_results'];
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
    disp('Your output folder did not exist, a new folder was created.')
end

%% load data table or create new one 
% load existing data
file_name         = '\wp10_results_post_nav.mat';
file_path         = fullfile(result_folder, file_name);
if isfile(file_path)
    load(file_path)
end

% initialize data if non-existing
if ~exist('pt','var')
    pt=[];
    p=1; % initialize participant index   
else
    p=0; % dummy value 
end 
 
%% Block 2: Data analysis
for id=participant_start:participant_end
    %% set participant index
    if p~=1 
        % check if ID exists in data 
        p_ind = find([pt.id]==id); 
        if isempty(p_ind) % if not: append participant data    
            [~,n]=size(pt);
            p=n+1; clear n p_ind; 
        else % if yes: overwrite participant data
            p=p_ind; clear p_ind;
            fprintf('Data for participant %d already existed and is overwritten now.\n', id);
        end
    end 
    
    %% set individual input folder
    % input folder
    input_folder=[data_folder '\' num2str(id) '\S00' num2str(s)];
    if ~exist(input_folder, 'dir')
        fprintf('Folder for participant %d, session %d does not exist. Continue with next iteration.\n', id, s);
        continue;
    end

    %% read-in trial results file
    opts=detectImportOptions([input_folder,'\trial_results.csv']);
    opts=setvaropts(opts,'timestamp','TreatAsMissing','na','Type','datetime','InputFormat','MM/dd/uuuu hh:mm:ss aa');
    trial_data=readtable([input_folder,'\trial_results.csv'], opts); clear opts;
      
    %% set participant, group and session duration information 
    pt(p).id=id; 
    [pt(p).group,pt(p).sex]=setGroupSexInfo(pt(p).id);
    pt(p).session_duration=round(minutes(trial_data.timestamp(4,1) - trial_data.timestamp(1,1))); 

    % loop over trials
    for k=1:size(trial_data,1)
        %% general info
        % trial
        pt(p).trial(k).trial=trial_data.trial_num(k,1); 
        % correct items 
        pt(p).trial(k).lm_MB=trial_data.lm_MB{1,1};
        pt(p).trial(k).lm_MD=trial_data.lm_MD{1,1};
        pt(p).trial(k).lm_MF=trial_data.lm_MF{1,1};
        pt(p).trial(k).lm_MH=trial_data.lm_MH{1,1};
        pt(p).trial(k).lm_MJ=trial_data.lm_MJ{1,1}; 
        pt(p).trial(k).obj_MA=trial_data.obj_MA{1,1};
        pt(p).trial(k).obj_MC=trial_data.obj_MC{1,1};
        pt(p).trial(k).obj_MI=trial_data.obj_MI{1,1};
        % dummy for chosen items 
        pt(p).trial(k).obj_1="999";
        pt(p).trial(k).obj_2="999";
        pt(p).trial(k).obj_3="999";
        pt(p).trial(k).obj_4="999";
        pt(p).trial(k).obj_5="999";
        
        %% time analysis
        b=trial_data.end_time(k,1); a=trial_data.start_time(k,1);
        pt(p).trial(k).time=computeTime(a,b); clear a b; 
        % fprintf('Time analysis done for %d, file no %d.\n', id, k);
        
        %% performance analysis
        if k==1 % SHAPE RECOGNITION 
            % score
            if trial_data.suc_1{k,1}=="True"
                points=1;
            elseif trial_data.suc_1{k,1}=="na"
                points=999;
            else
                points=0;
            end
            pt(p).trial(k).score=points;
            
            % save item info
            if trial_data.suc_1{k,1}=="na"
                pt(p).trial(k).obj_1="999"; 
            else
                pt(p).trial(k).obj_1=trial_data.obj_1{k,1};
            end
            % fprintf('Shape recognition done for %d, file no %d.\n', id, k);
            
        elseif k==2 % LANDMARK RECOGNITION      
            % get chosen items 
            obj_1=trial_data.obj_1{k,1}(4:end); obj_1=strsplit(obj_1, '_'); 
            obj_2=trial_data.obj_2{k,1}(4:end); obj_2=strsplit(obj_2, '_');
            obj_3=trial_data.obj_3{k,1}(4:end); obj_3=strsplit(obj_3, '_');
            obj_4=trial_data.obj_4{k,1}(4:end); obj_4=strsplit(obj_4, '_');
            obj_5=trial_data.obj_5{k,1}(4:end); obj_5=strsplit(obj_5, '_');
            obj_list=[obj_1; obj_2; obj_3; obj_4; obj_5];
            
            % score
            pt(p).trial(k).score=sum(contains(obj_list(:,2),'corr'))/size(obj_list,1);
            
            % save item info 
            pt(p).trial(k).obj_1=trial_data.obj_1{k,1};
            pt(p).trial(k).obj_2=trial_data.obj_2{k,1};
            pt(p).trial(k).obj_3=trial_data.obj_3{k,1};
            pt(p).trial(k).obj_4=trial_data.obj_4{k,1};
            pt(p).trial(k).obj_5=trial_data.obj_5{k,1}; 
            % fprintf('Landmark recognition done for %d, file no %d.\n', id, k);
            
        elseif k==3 % GOAL RECOGNITION
            % get chosen items
            corr_list={ trial_data.obj_MA{k,1}, trial_data.obj_MC{k,1}, trial_data.obj_MI{k,1} };
            obj_1=trial_data.obj_1{k,1};   
            obj_2=trial_data.obj_2{k,1};
            obj_3=trial_data.obj_3{k,1};
            obj_list={ obj_1; obj_2; obj_3 };
            
            % score 
            pt(p).trial(k).score=sum(contains(obj_list, corr_list))/size(obj_list,1);

            % save item information
            pt(p).trial(k).obj_1=obj_1;
            pt(p).trial(k).obj_2=obj_2;
            pt(p).trial(k).obj_3=obj_3;
            % fprintf('Goal recognition done for %d, file no %d.\n', id, k);
           
        else % POSITIONING: done externally with GMDA software
            pt(p).trial(k).score=999;
        end

        clear points obj* corr*; 
    end
    
    save(file_path, 'pt');
    fprintf('Analysis done for %d.\n', id);
    clear k trial_data input_folder; 
    p=0; % dummy value
end

%% Block 3: Write data to xlsx file
% [data_folder]  = setInputPath();
% writePostTableToXLSX(data_folder); 

clear; 
clear; close all; clc; format long;
addpath(genpath(pwd)) % add subfolder functions to path

%% Starmaze Data Processing
% @ date 2019-10-01 @ author Deetje Iggena
% @ date 2020-11-06 updated by @ Patrizia Maier & now tracked by git
% for Starmaze version WP10 Frankfurt
% requires Matlab 2021a or later (for shortestpath())
% requires Signal Processing Toolbox for downsample() and dtw()

% The script requires .csv files as input.
% One tracker file per trial with timestamp, x- and y-coordinates for movememt
% and z-coordinates for rotation.
% One trial_results file with general information (id, session, condition).

% Block 1: Set up input/output folders, Starmaze and Practise environment
% Block 2: Data preparation 
% Block 3: Data analysis
% Block 4: Write data to xlsx file

%% Block 1: Set up input and output folders, starmaze and practise environment 
%% data input folder and participant information
[data_folder]   = setInputPath(); % provide folder with all raw data
[participant_start,participant_end] = setParticipants(); % provide participant range
n_sessions      = 3; 

%% result folder
result_folder=[data_folder '\WP10_results'];
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
    disp('Your output folder did not exist, a new folder was created.')
end

%% load data table or create new one 
% load existing data
file_name       = '\wp10_results_navigation.mat';
file_path       = fullfile(result_folder,file_name);
if isfile(file_path)
    load(file_path)
end

% initialize if non-existing
if ~exist('sm','var')
    sm=[]; 

    %% set up Starmaze environment 
    % coordinates of min/max values
    values=table2array(readtable('wp10_values.csv'));
    [sm.coord.xmin, sm.coord.xmax, sm.coord.ymin,sm.coord.ymax]=setMinMaxValues(values);

    % coordinates of start positions (normalized!)
    start=table2array(readtable('wp10_start.csv'));
    [sm.coord.start_x,sm.coord.start_y]=setStartValues(start,sm.coord.xmin,sm.coord.xmax,sm.coord.ymin,sm.coord.ymax);

    % coordinates of goal positions (normalized!)
    goal=table2array(readtable('wp10_goal.csv'));
    [sm.coord.goal_x,sm.coord.goal_y]=setGoalValues(goal,sm.coord.xmin,sm.coord.xmax,sm.coord.ymin,sm.coord.ymax);
    [sm.coord.goal_x_in_alleys,sm.coord.goal_y_in_alleys]=setGoalMatrix(sm.coord.goal_x,sm.coord.goal_y);
    
    % coordinates of alley corners (normalized!)
    alley_x=table2array(readtable('wp10_alley_x.csv'));
    [n_corners,n_alleys]=size(alley_x);
    for i_alley=1:n_alleys
        for i_corner=1:n_corners
            alley_x(i_corner,i_alley)=setNormalizedValues(alley_x(i_corner,i_alley),sm.coord.xmin,sm.coord.xmax);
        end
    end
    alley_y=table2array(readtable('wp10_alley_y.csv'));
    for i_alley=1:n_alleys
        for i_corner=1:n_corners
            alley_y(i_corner,i_alley)=setNormalizedValues(alley_y(i_corner,i_alley),sm.coord.ymin,sm.coord.ymax);
        end
    end

    % coordinates of combined pentagon (normalized!) and central polyshape
    pentagon_x=table2array(readtable('wp10_pentagon_x.csv'));
    pentagon_y=table2array(readtable('wp10_pentagon_y.csv'));
    [sm.coord.central_x,sm.coord.central_y,sm.coord.central_poly,pentagon_x,pentagon_y]=setPentagonValues(alley_x,alley_y,pentagon_x,pentagon_y,...
        sm.coord.xmin,sm.coord.xmax,sm.coord.ymin,sm.coord.ymax);

    % create other polyshapes
    % alley polyshape
    [sm.coord.alley_full_x,sm.coord.alley_full_y,sm.coord.alley_poly,...
        sm.coord.alley_half_out_x,sm.coord.alley_half_out_y,sm.coord.alley_poly_out,...
        sm.coord.alley_half_in_x,sm.coord.alley_half_in_y,sm.coord.alley_poly_in]=setAlleyPolyshape(alley_x,alley_y);
    % rectangle polyshape
    [sm.coord.rec_x,sm.coord.rec_y,sm.coord.rec_poly]=setRectPolyshape(n_alleys,alley_x,alley_y,pentagon_x,pentagon_y);
    % triangle polyshape
    [sm.coord.tri_x,sm.coord.tri_y,sm.coord.tri_poly]=setTriPolyshape(n_alleys,alley_x,alley_y,pentagon_x,pentagon_y);
    % joint polyshape
    sm.coord.full_poly=[sm.coord.alley_poly_out{1,1} sm.coord.alley_poly_in{1,1}...
        sm.coord.alley_poly_out{1,2} sm.coord.alley_poly_in{1,2}...
        sm.coord.alley_poly_out{1,3} sm.coord.alley_poly_in{1,3}...
        sm.coord.alley_poly_out{1,4} sm.coord.alley_poly_in{1,4}...
        sm.coord.alley_poly_out{1,5} sm.coord.alley_poly_in{1,5} sm.coord.central_poly];
    
    % compute random x-/y-coordinate distribution in Starmaze
    [sm.coord.random_x, sm.coord.random_y]=computeRandomLocations(sm.coord.full_poly, 1000); 
    
    % compute final distance between locations and random x-/y-coordinates   
    [sm.coord.final_distance_distribution]=computeDistanceDistribution(sm.coord.goal_x_in_alleys,...
        sm.coord.goal_y_in_alleys, sm.coord.random_x, sm.coord.random_y); 
    
    % create graph
    % for automated shortest path calculation (requires Matlab 2021a)
    [sm.coord.graph,sm.coord.graph_x,sm.coord.graph_y]=setGraph(sm.coord.start_x,sm.coord.start_y,...
        sm.coord.tri_x,sm.coord.tri_y,sm.coord.goal_x,sm.coord.goal_y);

    % add (ordered) information
    sm.coord.goal_names=["MA", "MC", "MI"];
    sm.coord.start_names=["A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "X"];
    sm.coord.alley_names=["A" "B" "C" "D" "E" "F" "G" "H" "I" "J"];

    % create test figure plot
    % plotTestFigure("s_maze",sm.coord.full_poly,sm.coord.goal_x,sm.coord.goal_y,sm.coord.start_x,sm.coord.start_y,sm.coord.goal_names,sm.coord.graph,sm.coord.graph_x,sm.coord.graph_y);

    clear i j distribution alley_x alley_y pentagon_x pentagon_y i_alley n_alleys i_corner n_corners start goal values; 
    
    %% set up Practise environment (for motor control task)
    % coordinates of min/max values
    practise_values=table2array(readtable('wp10_practise_values.csv'));
    [sm.coord.practise.xmin,sm.coord.practise.xmax,sm.coord.practise.ymin,sm.coord.practise.ymax]=setMinMaxValues(practise_values);

    % coordinates of start position (normalized!)
    practise_start=table2array(readtable('wp10_practise_start.csv'));
    [sm.coord.practise.start_x,sm.coord.practise.start_y]=setStartValues(practise_start,sm.coord.practise.xmin,sm.coord.practise.xmax,sm.coord.practise.ymin,sm.coord.practise.ymax);

    % coordinates of goal positions (normalized!)
    practise_goal=table2array(readtable('wp10_practise_goal.csv'));
    [sm.coord.practise.goal_x,sm.coord.practise.goal_y]=setGoalValues(practise_goal,sm.coord.practise.xmin,sm.coord.practise.xmax,sm.coord.practise.ymin,sm.coord.practise.ymax);

    % coordinates of alley corners (normalized!)
    practise_alley_x=table2array(readtable('wp10_practise_x.csv'));
    [n_corners,~] = size(practise_alley_x);
    for i_corner=1:n_corners
        practise_alley_x(i_corner,1)=setNormalizedValues(practise_alley_x(i_corner,1),sm.coord.practise.xmin,sm.coord.practise.xmax);
    end
    practise_alley_y=table2array(readtable('wp10_practise_y.csv'));
    for i_corner=1:n_corners
        practise_alley_y(i_corner,1)=setNormalizedValues(practise_alley_y(i_corner,1),sm.coord.practise.ymin,sm.coord.practise.ymax);
    end

    % create polyshape
    sm.coord.practise.practise_poly=polyshape(practise_alley_x(:),practise_alley_y(:));

    % add (ordered) information
    sm.coord.practise.practise_goal_names=["1", "2", "3", "4", "5", "6", "7", "8" , "9", "10" ];
    sm.coord.practise.practise_start_names="Player_P0";

    % create test figure plot
    % plotTestFigure("p_maze",sm.coord.practise.practise_poly,sm.coord.practise.goal_x,sm.coord.practise.goal_y,sm.coord.practise.start_x,sm.coord.practise.start_y,sm.coord.practise.practise_goal_names,sm.coord.graph,sm.coord.graph_x,sm.coord.graph_y);

    clear practise_alley_x practise_alley_y i_corner n_corners practise_start practise_goal practise_values; 
    
    %% initialize participant index
    p=1; 
else
    p=0; % default value 
end

%% Block 2: Data preprocessing
for id=participant_start:participant_end
tic;     
    % set participant index
    if p~=1 
        % check if ID exists in data 
        p_ind = find([sm.participant.id]==id); 
        if isempty(p_ind) % if not: append participant data    
            [~,n]=size(sm.participant);
            p=n+1; clear n p_ind; 
        else % if yes: overwrite participant data
            p=p_ind; clear p_ind;
            fprintf('Data for participant %d already existed and is overwritten now.\n', id);
        end
    end 
    
    % loop over sessions 
    for s=1:n_sessions        
        %% set individual input and output folder
        % input folder
        input_folder=[data_folder '\' num2str(id) '\S00' num2str(s)]; 
        if ~exist(input_folder, 'dir')
            fprintf('Folder for participant %d, session %d does not exist. Continue with next iteration.\n', id, s);
            continue;
        end
        
        % output folder (only for trial plots)
        output_folder=[data_folder '\' num2str(id) '\plots'];
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
            fprintf('Your output folder for participant %d did not exist, a new folder was created.\n', id);
        end
        
        %% read-in trial results file
        opts=detectImportOptions([input_folder, '\trial_results.csv']);
        opts=setvaropts(opts,'timestamp','InputFormat','MM/dd/uuuu hh:mm:ss aa');
        trial_data=readtable([input_folder, '\trial_results.csv'], opts); clear opts; 
        
        % store participant information
        sm.participant(p).id=id;
        [sm.participant(p).group,sm.participant(p).sex]=setGroupSexInfo(sm.participant(p).id);
        sm.participant(p).session(s).session=trial_data.session_num(1,1);
        
        %% read-in log file (only once in session 1)
        if s==1
            opts=detectImportOptions([input_folder, '\log.csv'], 'VariableNamesLine', 1, 'Delimiter', ',');
            opts.DataLines=[2,Inf]; 
            opts.SelectedVariableNames={'message'};
            log_data=table2cell(readtable([input_folder, '\log.csv'], opts));
            log_data=log_data(contains(log_data,'ID is'));
            [sm.participant(p).rand_dict]=setRandomizationDict(log_data, id);
            clear log_data opts; 
        end 
        
        %% get individual trial tracker file names 
        d=dir(fullfile(input_folder, 'trackedpoint_*.csv'));
        files={d.name}; clear d;
        
        % loop over trials 
        for k=1:numel(files)
            name=files{k};
            
            % for session 3 only process trial 2 (motor control task)
            if s==3
                pattern=("_T001"|"_T003"|"_T004"|"_T005"|"_T006"|"_T007");
                if contains(name, pattern)
                    continue
                end
            end
            
            % read-in and clean trial tracker data
            data=readtable(fullfile(input_folder, name));
            data=cleanTrialTrackerFile(data,s); 
                    
            % extract data 
            t=data.time; % time
            x=data.pos_x; % x coordinates
            y=data.pos_z; % y coordinates
            r=data.rot_y; R=deg2rad(r); r=unwrap(R); % yaw rotation
            if s==3
                events=data.trialEvent; 
            end
            clear data R;  
            
            % spatial normalization
            if s==3 % practise maze
                x=setNormalizedValues(x,sm.coord.practise.xmin,sm.coord.practise.xmax); 
                y=setNormalizedValues(y,sm.coord.practise.ymin,sm.coord.practise.ymax);
            else % star maze
                x=setNormalizedValues(x,sm.coord.xmin,sm.coord.xmax); 
                y=setNormalizedValues(y,sm.coord.ymin,sm.coord.ymax);
            end
            % save end points
            sm.participant(p).session(s).trial(k).x_n=x(end); sm.participant(p).session(s).trial(k).y_n=y(end);
            
            % temporal normalization
            % set new sampling rate (default 0.05 corresponds to 20 frames per second)          
            new_sampling_rate=0.05; 
            [ti,xi,yi,ri]=temporalNormalization(new_sampling_rate,t,x,y,r); 
            clear t x y r new_sampling_rate; 
                                  
            %% get single trial info from trial_results
            sm.participant(p).session(s).session_duration=round(minutes(trial_data.timestamp(numel(files),1) - trial_data.timestamp(1,1))); 
            
            sm.participant(p).session(s).trial(k).block=trial_data.block_num(k,1);
            sm.participant(p).session(s).trial(k).trial=trial_data.trial_num(k,1);
            sm.participant(p).session(s).trial(k).trial_in_block=trial_data.trial_num_in_block(k,1);
            
            sm.participant(p).session(s).trial(k).feedback=string(trial_data.trial_feedback(k,1));
               
            condition=string(trial_data.trial_type(k,1));
            sm.participant(p).session(s).trial(k).condition=setTrialCondition(condition,sm.participant(p).session(s).trial(k).feedback);
            clear condition; 
            
            sm.participant(p).session(s).trial(k).goal_identity=string(trial_data.trial_goal_identity(k,1));
            
            [sm.participant(p).session(s).trial(k).goal_x,sm.participant(p).session(s).trial(k).goal_y,...
                sm.participant(p).session(s).trial(k).goal_i,~,...
                sm.participant(p).session(s).trial(k).goal_alley]=setTrialGoalLocation(char(string(trial_data.trial_goal(k,1))),...
                sm.coord.goal_x,sm.coord.goal_y,sm.coord.goal_names,sm.coord.alley_names);
            
            start_name=char(trial_data.trial_player(k,1)); start_s=string(start_name(end)); 
            [sm.participant(p).session(s).trial(k).start_i]=setTrialStartLocation(start_s,sm.coord.start_names);
            clear start_name start_s; 
            
            %% For all normal navigation trials (i.e., not motor control task)
            if sm.participant(p).session(s).trial(k).condition~="practise"
                %% compute support variables depending on this trial's settings
                % ideal path coordinates & ideal path length
                % Requires Matlab 2021a for 'shortestpath' function and
                % 'interparc' function by John D'Errico (Matlab File Exchanger)
                [x_line, y_line, x_line_chosen, y_line_chosen, ~, ~,...
                    ~, ~, ~, ~, ~, ideal_path_length_to_target, ideal_path_length_to_chosen, ~]=computeStartDependentIdealVariables(...
                    sm.coord.graph, sm.coord.graph_x, sm.coord.graph_y,...
                    sm.participant(p).session(s).trial(k).start_i, sm.participant(p).session(s).trial(k).goal_i,...
                    sm.participant(p).session(s).trial(k).x_n, sm.participant(p).session(s).trial(k).y_n,...
                    sm.coord.alley_full_x, sm.coord.alley_full_y, sm.coord.rec_x, sm.coord.rec_y, ...
                    sm.coord.central_poly, sm.coord.full_poly);
                               
                % interpolate ideal path data for further analysis
                % using 'interparc' function by John D'Errico (Matlab File Exchanger)
                [xi_al,yi_al]=interpolateData(x_line,y_line,ideal_path_length_to_target);
                                
%                 % test plot
%                 figure; plot(sm.coord.full_poly); hold on; 
%                 plot(xi, yi, 'k-', x_line, y_line, 'k+', xi_al, yi_al, 'b-',...
%                     x_line_chosen, y_line_chosen, 'm--');
%                 xlim([0 1]); ylim([0 1]); hold off; 
                               
                %% Block 3: Data analysis, i.e. calculcation of variables     
                %% accuracy analysis (only for probe trials)
                % compute chosen goal location
                [~, ~, sm.participant(p).session(s).trial(k).chosen_alley_i, ~]=computeChosenGoals(...
                    sm.participant(p).rand_dict,char(trial_data.chosen_goal(k,1)),...
                    sm.coord.alley_poly,sm.coord.rec_poly,sm.coord.tri_poly,sm.coord.alley_names,sm.coord.goal_names,...
                    sm.participant(p).session(s).trial(k).x_n,sm.participant(p).session(s).trial(k).y_n);

                % evaluate chosen goal location
                if sm.participant(p).session(s).trial(k).feedback=="false"  
                    % FINAL DISTANCE to TARGET 
                    final_distance=computeDistance(sm.participant(p).session(s).trial(k).goal_x,...
                        sm.participant(p).session(s).trial(k).x_n,sm.participant(p).session(s).trial(k).goal_y,...
                        sm.participant(p).session(s).trial(k).y_n);
                   
                    % MEMORY SCORE (final distance in relation to distribution from random points)
                    sm.participant(p).session(s).trial(k).memory_score=computeMemoryScore(sm.coord.final_distance_distribution,...
                        final_distance,sm.participant(p).session(s).trial(k).goal_i,sm.participant(p).session(s).trial(k).goal_alley); 
                    
                    % CORRECT FINAL ALLEY 
                    sm.participant(p).session(s).trial(k).correct_final_alley=...
                        sm.participant(p).session(s).trial(k).goal_alley==sm.participant(p).session(s).trial(k).chosen_alley_i;
                else
                    % default values
                    sm.participant(p).session(s).trial(k).memory_score=999; 
                    sm.participant(p).session(s).trial(k).correct_final_alley=999; 
                end
                % fprintf('Accuracy analysis done for %d, session %d, file no %d.\n', id, s, k);
                
                %% time analysis
                % TIME
                sm.participant(p).session(s).trial(k).time=computeTime(ti(1),ti(end));  
                                             
                %% standard coordinate analysis using x-/y-coordinates
                % PATH LENGTH 
                path_length=computePathLength(xi,yi); 
                
                % EXCESS PATH LENGTH (to chosen target) 
                sm.participant(p).session(s).trial(k).excess_path_length=path_length - ideal_path_length_to_chosen;
                                 
                % AVERAGE DISTANCE to GOAL / PROXIMITY
                [distance_to_goal, ~]=computeTargetProximity(xi,yi,...
                    sm.participant(p).session(s).trial(k).goal_x,sm.participant(p).session(s).trial(k).goal_y); 
                
                % EXCESS AVERAGE DISTANCE to GOAL / PROXIMITY 
                % approximation for 'ideal value' on shortest path with equally spaced steps
                [ideal_distance_to_goal, ~]=computeTargetProximity(xi_al,yi_al,...
                    sm.participant(p).session(s).trial(k).goal_x,sm.participant(p).session(s).trial(k).goal_y);               
                % excess difference 
                % value < 0 is possible and indicates more time was spend close to goal than necessary
                sm.participant(p).session(s).trial(k).excess_distance_to_goal=distance_to_goal - ideal_distance_to_goal; 
                
                % fprintf('Standard time, path & distance analysis done for %d, session %d, file no %d.\n', id, s, k);
          
                %% rotation analysis using z-rotation
                % TOTAL ROTATION
                % cumulative absolute change in yaw rotation (r)
                 [sm.participant(p).session(s).trial(k).total_rotation]=computeRotation(ri); 
                
                % INITIAL ROTATION in this trial's START AREA (alley plus triangle)
                % get rotation index (different method for inner/outer starts)
                if mod(sm.participant(p).session(s).trial(k).start_i,2)
                    [~, ~, ~, rot_index]=computePresence(sm.participant(p).session(s).trial(k).start_i,xi,yi,...
                        sm.coord.alley_poly, sm.coord.tri_poly, sm.participant(p).session(s).trial(k).time);
                else 
                    [~, ~, ~, rot_index]=computePresencePentagon(sm.participant(p).session(s).trial(k).start_i,xi,yi,...
                        sm.coord.rec_poly, sm.coord.tri_poly, sm.participant(p).session(s).trial(k).time); 
                end 
                [rot_index]=computeFirstSegment(rot_index); 
                % compute initial rotation
                % cumulative absolute change in yaw rotation (r)
                [sm.participant(p).session(s).trial(k).initial_rotation]=computeRotation(ri(rot_index));
                % initial_rotation_velocity=sm.participant(p).session(s).trial(k).initial_rotation/(numel(ri(rot_index))-1);            
                clear rot_index;
                
                % fprintf('Rotation analysis done for %d, session %d, file no %d.\n', id, s, k);   
                  
                %% set marker for excluded trials
                % criteria: timeout, or no movement, or very short trial time
                sm.participant(p).session(s).trial(k).exclude_trial_matlab=0;
                if sm.participant(p).session(s).trial(k).chosen_alley_i==999
                    sm.participant(p).session(s).trial(k).exclude_trial_matlab=1;
                    fprintf('Trial %d marked for exclusion due to timeout.\n',k);
                elseif (path_length<=0.1 || sm.participant(p).session(s).trial(k).total_rotation==0 || ...
                        sm.participant(p).session(s).trial(k).time<3)
                    sm.participant(p).session(s).trial(k).exclude_trial_matlab=1;
                    fprintf('Trial %d marked for exclusion due lack of movement or trial time < 3 sec.\n',k);
                end
                
                %% create plot  
                plotTrialTrack(sm.participant(p).id, s,...
                    sm.participant(p).session(s).trial(k).trial,...
                    sm.participant(p).session(s).trial(k).condition,...
                    sm.participant(p).session(s).trial(k).memory_score,...
                    sm.participant(p).session(s).trial(k).time,...
                    sm.participant(p).session(s).trial(k).excess_path_length,...
                    sm.participant(p).session(s).trial(k).excess_distance_to_goal,...
                    sm.participant(p).session(s).trial(k).initial_rotation,...
                    sm.coord.full_poly,xi,yi,x_line,y_line,x_line_chosen,y_line_chosen,...
                    sm.participant(p).session(s).trial(k).goal_x,sm.participant(p).session(s).trial(k).goal_y,...
                    output_folder);

            else 
                %% For motor control task               
                %% compute support variables depending on the trial's settings
                % ideal path coordinates & length
                x_line_motor=[sm.coord.practise.start_x; sm.coord.practise.goal_x]; 
                y_line_motor=[sm.coord.practise.start_y; sm.coord.practise.goal_y];
                ideal_motor_path_length=computePathLength(x_line_motor,y_line_motor);
                               
                %% time analysis
                % TIME
                sm.participant(p).session(s).trial(k).time=computeTime(ti(1),ti(end));

                %% standard coordinate analysis using x-/y-coordinates and z-rotation
                % PATH LENGTH to all targets
                path_length=computePathLength(xi,yi); 
                              
                % EXCESS PATH LENGTH to all targets
                sm.participant(p).session(s).trial(k).excess_path_length=path_length-ideal_motor_path_length;

                % TOTAL ROTATION
                [sm.participant(p).session(s).trial(k).total_rotation]=computeRotation(ri);

                % fprintf('Motor control analysis done for %d, session %d, file no %d.\n', id, s, k);
                
                %% set marker for excluded trials to zero
                sm.participant(p).session(s).trial(k).exclude_trial_matlab=0;
                
                %% create plot  
                plotMotorControlTrack(sm.participant(p).id, s, sm.participant(p).session(s).trial(k).trial,...
                    sm.participant(p).session(s).trial(k).time,...
                    sm.participant(p).session(s).trial(k).excess_path_length,...
                    sm.participant(p).session(s).trial(k).total_rotation,...
                    sm.coord.practise.practise_poly,sm.coord.practise.goal_x,sm.coord.practise.goal_y,...
                    sm.coord.practise.start_x,sm.coord.practise.start_y,...
                    sm.coord.practise.practise_goal_names,xi,yi,x_line_motor,y_line_motor,output_folder);
                clear x_line_motor y_line_motor;     
            end
            
            clear x* y* ti ri ideal_path* origin* name;             
        end
        
        clear files trial_data input_folder; 
        fprintf('Analysis done for %d, session %d.\n', id, s);
    end
    
    p=0; % default value
    save(file_path, 'sm'); 
    toc;
end

%% Block 4: Write data to xlsx file
% [data_folder]  = setInputPath();
% writeNavTableToXLSX(data_folder); 

clear; 
function data=cleanTrialTrackerFile(data,session)
% cleanTrialTrackerFile Data cleaning of trial tracker files. 
%
% Input: 
% data is input trial tracker data (table).
% session is session number (integer).
%
% Returns: data as cleaned trial tracker data (table). 

% clean rows
data(1:2,:)=[]; % remove first two rows to get rid of info from last trial
data(strcmp(data.gameIsPaused,'True'),:)=[]; % remove rows when gameIsPaused
if session~=3 % skip this for motor control task
    data(data.trialEvent==1,:)=[]; % remove rows when trialEvent == 1 (i.e. after goal was found)
end
% clean columns 
data.pos_y=[]; 
data.rot_x=[]; 
data.rot_z=[]; 
data.gameIsPaused=[]; 

end
    

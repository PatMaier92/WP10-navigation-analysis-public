function [start_int]=setTrialStartLocation(start_string, start_names)
% setTrialStartLocation: Returns start position information for this trial 
% for Starmaze WP10. 
%
% Input: 
% start_string ist start position (string).
% start_names is ordered string array of starting points. 
%
% Returns: start position information (integer).

% integer
start_int=find(contains(start_names, start_string)); 

% correct for ego: same start_int as in learning 
if start_int==11 % X
    start_int=7; % G
end

% correct for false input 
if isempty(start_int)
    disp('Unknown starmaze input information in trialStart.m. Set to 999.');
    start_int=999; 
end

end
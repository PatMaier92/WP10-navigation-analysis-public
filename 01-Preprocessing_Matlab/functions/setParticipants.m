function [participant_a,participant_z]=setParticipants()
% setParticipants: Takes two integers as input. 
% Returns them as seperate variables. 
%
% Input: 
% e.g. 1111, 1120
%
% Returns: participant_a, participant_z.

% first subject
isValid = false;
while ~isValid
    participant_a = str2double(input('Enter first participant: ','s'));
    if isnan(participant_a)
        clear participant_a;
    else
        isValid = true;
    end
end

% last subject
isValid = false;
while ~isValid
    participant_z = str2double(input('Enter last participant: ','s'));
    if isnan(participant_z)
        clear participant_z; 
    else
        isValid = true;
    end
end

end


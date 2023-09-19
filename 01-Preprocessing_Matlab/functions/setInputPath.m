function [data_path]=setInputPath()
% setInputPath: Takes input path, checks for validity and stores it. 
%
% Input: 
% e.g D:\Temp Home Office\Analysis\Functions.
% 
% Returns: data_path is saved path. 

isValid = false;
while ~isValid
    data_path = input('Enter input folder, where all data is located: ','s');
    if ~exist(data_path, 'dir')
        disp('Your input folder does not exist. Enter a valid folder.')
    else
        isValid = true;
    end
end

end


function [rand_dict]=setRandomizationDict(log_data, subject)
% setRandomizationDict Preprocessing of log data with information on
% randomization of goals and starts.
%
% Input: cleaned log data (cell).
%
% Returns: rand_dict (structure) contains randomization info.

if size(log_data,1) ~= 4
    disp('Error: Log file is too short. Please check the file.\n');
else
    for i=1:4
        line=split(log_data(i),' ');
        line=line(~cellfun('isempty',line)); 
        
        % id and key
        id=str2double(line{3});
        if id ~= subject
            disp('Error: ID mismatch detected during log data randomization check. Please check the raw data.\n');
            %break
        end
        key = line{5}; 
        
        % rest of information 
        inner_dict={};
        line=line(7:end); 
        for n=1:2:14
            k=line{n}; v=line{n+1};           
            inner_dict.(k)=v;
        end
        
        % add T1 recall and T2 goal order (hard-coded in Unity)
        if str2double(inner_dict.T1)==1
            inner_dict.T1_Recall=3;
            inner_dict.T2_Recall=2; 
        elseif str2double(inner_dict.T1)==2
            inner_dict.T1_Recall=1;
            inner_dict.T2_Recall=3;
        elseif str2double(inner_dict.T1)==3
            inner_dict.T1_Recall=2;
            inner_dict.T2_Recall=1;
        elseif str2double(inner_dict.T1)==0
            inner_dict.T1_Recall=0;
            inner_dict.T2_Recall=0;
        end
        
        rand_dict.(key)=inner_dict;
    end
end

end
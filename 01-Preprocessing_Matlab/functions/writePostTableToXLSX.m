function writePostTableToXLSX(data_folder)
% writePostTableToXLSX Write data summary table in xlsx format for post-navigational memory test. 
%
% Input: path to data folder 
%
% Returns: data summary table as xlsx file 

% load data 
file_path=[data_folder '\WP10_results']; 
file_name='\wp10_results_post_nav.mat';
full_file=fullfile(file_path, file_name); 
if isfile(full_file)
    load(full_file, 'pt'); 
else 
    disp('Your input data does not exist. Please check your data folder.'); 
    return; 
end 
data=pt; clear pt; 

% process data (from structure to table) 
temp=[];
session=4; 
[~,p]=size(data);
[~,r]=size(data(1).trial); 
for i=1:p
    g_data=table(repmat(data(i).id,r,1), repmat(string(data(i).group),r,1), repmat(string(data(i).sex),r,1),...
        repmat(session,r,1), repmat(data(i).session_duration,r,1),...
        'VariableNames',{'id' 'group' 'sex' 'session' 'duration'});
    t_data=struct2table(data(i).trial);
    gt_data=[g_data t_data];
    % add to data table 
    temp=vertcat(temp,gt_data);
    clear *_data;
end

% sort order 
temp=sortrows(temp,{'id','trial'});

% write data 
format='yymmdd'; date=datestr(now, format); 
file_name2=['wp10_post_navigation_data_' date '.xlsx']; 
writetable(temp,fullfile(file_path,file_name2));

end
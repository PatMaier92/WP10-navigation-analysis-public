function [index]=computeFirstSegment(original_index)
% computeFirstSegment: Determines index of first segment of presence in zone. 
% 
% Input: 
% original_index with all segments of presence in zone (boolean)
% 
% Returns:
% index of first segment of presence in zone (boolean) 

% default
index=false(numel(original_index),1); 

if all(original_index==1)
    % participant did not leave initial zone: set index to true
    index=true(numel(original_index),1); 
else 
    % find index(es) of segments in this zone
    i_start=find(diff(original_index) > 0); % start index(es)
    i_end=find(diff(original_index) < 0); % end index(es)
    if ~isempty(i_start) || ~isempty(i_end)
        
        % correct index(es) for missing start or end
        if length(i_start) < length(i_end) % add start
            i_start=[1; i_start];
        elseif length(i_start) > length(i_end)  % add end
            i_end(end+1)=length(index);
        end
        
        % create index of first segment
        index=[true(i_end(1),1); false(numel(original_index)-i_end(1),1)];
    end
end 

end 

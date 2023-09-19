function [fd_distribution]=computeDistanceDistribution(...
    set_x1, set_y1, set_x2, set_y2)
% computeDistanceDistribution Calculates Euclidean distance distribution
% between two sets of x-/y-coordinates of locations. 
%
% Input: 
% set_x1, set_y1 are x-/y-coordinates of first set, e.g. goal locations (float)
% set_x2, set_y2 are x-/y-coordinates of second set, e.g. random locations (float)
%
% Returns:  
% fd_distribution is a matrix with Euclidean distance distributions with 
% the same size as the first set (float) 
    
fd_distribution={};
for i=1:size(set_x1,1)
    for j=1:size(set_x1,2)
        % skip impossible rotations (x & y are zero)
        if set_x1(i,j)==0
            distribution=[];
        else
            % set default values
            distribution=zeros(numel(set_x2),1);
            % compute final distance distribution
            for m=1:numel(set_x2)
                distribution(m)=computeDistance(set_x1(i,j),...
                    set_x2(m), set_y1(i,j), set_y2(m));
            end
        end
        % save
        fd_distribution{i,j}=sort(distribution, 'ascend');
    end
end

end
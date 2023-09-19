function [avg_distance, total_distance]=computeTargetProximity(x, y, x_target, y_target)
% computeTargetProximity Calculates target proximity as averaged 
% sum of distance traveled between points.
%
% Input: 
% x, y are vectors with x-/y-coordinate points, i.e. a path. 
% x_target, y_target are target x-/y-coordinates.
%
% Returns:  
% avg_distance (float) and total_distance (float).

total_distance=0; 
for i=1:numel(x)-1
    total_distance=total_distance+computeDistance(x(i),x_target,y(i),y_target); 
end

avg_distance=total_distance/(numel(x)-1);
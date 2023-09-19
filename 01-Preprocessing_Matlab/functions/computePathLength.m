function path_length=computePathLength(x_line, y_line)
% computePathLength Calculates sum of distance traveled between points.
%
% Input: 
% x_line, y_line are vectors with x-/y-coordinate points, i.e. a path. 
%
% Returns:  
% Path length as summation of distance between points (float).

path_length=0; 
for i=1:length(x_line)-1
    path_length=path_length+computeDistance(x_line(i),x_line(i+1), y_line(i),y_line(i+1)); 
end
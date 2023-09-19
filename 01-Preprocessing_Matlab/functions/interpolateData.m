function [xi, yi]=interpolateData(x_line, y_line, ideal_distance)
% interpolateData: Interpolating data for Starmaze WP10 depending
% on path line vectors. This method takes into account distance values 
% to create equally spaced interpolated values. 
% Using the 'interparc' function by John D'Errico (Matlab File Exchanger). 
%
% Input:
% x_line,y_line are vectors with x-/y-coordinates
% ideal_distance is ideal path distance values
%
% Returns:
% xi,yi are vectors with interpolated x-/y-coordinates

% interpolate
crit=round(ideal_distance*1000);
xy=interparc(crit,x_line,y_line,'linear');
xi=xy(:,1); yi=xy(:,2);

end
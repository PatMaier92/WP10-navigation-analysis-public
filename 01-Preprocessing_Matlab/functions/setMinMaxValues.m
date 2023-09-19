function [xmin,xmax,ymin,ymax]=setMinMaxValues(values)
% setMinMaxValues: Takes Starmaze WP10 minimum and maximum
% x-,y-coordinates and returns them as single coordinate variables. 
%
% Input: 
% values are an array which is read-in from a csv. file. 
%
% Returns: xmin,xmax,ymin,ymax.

xmin=values(1,1);
xmax=values(2,1);

ymin=values(1,2);
ymax=values(2,2);

end
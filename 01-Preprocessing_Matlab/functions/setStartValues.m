function [start_x,start_y]=setStartValues(start,xmin,xmax,ymin,ymax)
% setStartValues: Normalizes start positions in Starmaze WP10 based on min/max
% xy-coordinates. 
%
% Input: 
% start is xy-coordinates of start position.
% xmin,xmax,ymin,ymax are minimum, maximum xy-coordinates.
%
% Returns: start_x,start_y.  

[row,col]=size(start);
for r=1:row
    start_x(r,1)=setNormalizedValues(start(r,1),xmin,xmax);
    start_y(r,1)=setNormalizedValues(start(r,2),ymin,ymax);
end

end
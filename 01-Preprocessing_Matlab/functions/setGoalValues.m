function [goal_x,goal_y]=setGoalValues(goal,xmin,xmax,ymin,ymax)
% setGoalValues: Normalizes goal location positions in Starmaze WP10 based on min/max
% xy-coordinates. 
%
% Input: 
% goal is xy-coordinates of goal location.
% xmin,xmax,ymin,ymax are minimum, maximum xy-coordinates.
%
% Returns: goal_x,goal_y.

[row,col]=size(goal);
for r=1:row
    goal_x(r,1)=setNormalizedValues(goal(r,1),xmin,xmax);
    goal_y(r,1)=setNormalizedValues(goal(r,2),ymin,ymax);
end

end
function [goal_x_in_alleys,goal_y_in_alleys]=setGoalMatrix(goal_x,goal_y)
% setGoalMatrix: Create matrix with hypothetical goal locations in each
% alley (by rotating goal locations)
% Input: 
% goal is xy-coordinates of all goal locations.
% polygon is polygon figure for plotting. 
%
% Returns: goal_x_in_alleys, goal_y_in_alleys.

goal_x_in_alleys=zeros(3,5);
goal_y_in_alleys=zeros(3,5);
x_center = 0.5; y_center = 0.5;

% for each rotation theta
for i=0:4
    % set rotation matrix
    theta = -360/5*i;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    
    % for each goal
    for j=1:3
        % set goal alley integer
        if j==1
            alley=1;
        elseif j==2
            alley=3;
        elseif j==3
            alley=9;
        end
        
        % rotate
        v = [goal_x(j) ; goal_y(j)];
        center = repmat([x_center; y_center], 1, length(v));
        vo = R*(v - center) + center;
        
        % store
        integer = (alley + 1)/2 + i;
        if integer > 5
            integer = integer - 5;
        end
        goal_x_in_alleys(j,integer) = vo(1,1);
        goal_y_in_alleys(j,integer) = vo(2,1);
    end
end
    
%     % test plot
%     figure; plot(sm.coord.full_poly); hold on;
%     plot(goal_x_in_alleys(1,:), goal_y_in_alleys(1,:), 'rx', ...
%         goal_x_in_alleys(2,:), goal_y_in_alleys(2,:), 'bx', ...
%         goal_x_in_alleys(3,:), goal_y_in_alleys(3,:), 'gx', ...
%         x_center, y_center, 'bo');
%     xlim([0 1]); ylim([0 1]); hold off;

end
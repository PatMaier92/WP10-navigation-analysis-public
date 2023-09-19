function plotMotorControlTrack(id, session, trial, time, excess_path, rotation, ...
    pract_polyshape, goal_x, goal_y, start_x, start_y,...
    pract_goal_locs, x, y, x_line, y_line, folder)
% plotMotorControlTrack: Creates plot for motor control trial.
%
% Input: Information for creating and naming the plot.
%
% Returns: Saves one nice plot.

ID=num2str(id); 
Session=num2str(session);
T=int2str(trial);
TI=num2str(round(time,1)); 
EP=num2str(round(excess_path,2)); 
RP=num2str(round(rotation,2)); 

wfig=figure('Position',[500 200 580 500],'visible','off');
plot(pract_polyshape);
axis([0 1 0 1]); xticks(0:0.1:1); yticks(0:0.1:1); 
hold on;
title({[ID ', Session: ' Session' ', Trial: ' T ', Condition: Motor Control']; ...
    ['time: ' TI, ', ex. path: ' EP, ', rotation: ' RP]});
for g=1:length(goal_x)
    viscircles([goal_x(g) goal_y(g)], 0.02);
end
viscircles([start_x start_y], 0.03);

% labels
for i=1:length(pract_goal_locs)
    text(goal_x(i)+0.02, goal_y(i), pract_goal_locs(i));
end 
text(start_x-0.05, start_y-0.05, 'Start');

% lines and legend 
line1=plot(x, y,'k -', 'LineWidth', 1);
line2=plot(x_line, y_line,'r -.', 'LineWidth', 1);
legend([line1 line2],{'actual path','ideal path'},'Location','north'); 
hold off;

% save plot
file_name =['Motor_Plot_' ID '_' Session '_' T '.jpeg'];
saveas(wfig, fullfile(folder, file_name)); 

end
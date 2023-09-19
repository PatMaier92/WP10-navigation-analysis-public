function plotTrialTrack(id, session, trial, condition,...
    memory_score, time, excess_path, excess_distance, rotation, ...
    polyshape, x, y, x_line, y_line, x_line_chosen, y_line_chosen,...
    goal_x, goal_y, folder)
% plotTrialTrack Creates track plots for each individual trial.
%
% Input: Information for creating and naming the plot.
%
% Returns: A nice trial track plot.

ID=num2str(id);
Session=num2str(session);
Trial=int2str(trial);
MS=num2str(round(memory_score,2)); 
TI=num2str(round(time,1)); 
EP=num2str(round(excess_path,2)); 
ED=num2str(round(excess_distance,2)); 
RP=num2str(round(rotation,2)); 

wfig=figure('visible','off');
plot(polyshape,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.1);
axis([0 1 0 1]); xticks(0:0.1:1); yticks(0:0.1:1); 
hold on;

viscircles([goal_x goal_y],0.01);

line1=plot(x,y,'k -', 'LineWidth', 1);
line2=plot(x_line,y_line,'r -.', 'LineWidth', 1);
line3=plot(x_line_chosen,y_line_chosen,'b .:', 'LineWidth', 1);
legend([line1 line2 line3],{'actual path','ideal path target','ideal path chosen'});

if condition=="main_learn"
    Type = 'Training'; 
elseif condition=="allo_ret"
    Type = 'Allo Probe';
elseif condition=="ego_ret"
    Type = 'Ego Probe';
elseif condition=="main_ret"
    Type = 'Training Probe';
else
    Type = ' (XXXXX)'; 
end
title({[ID ', Session: ' Session ', Trial: ' Trial ', Condition: ' Type]; ...
    ['memory: ', MS, ', time: ' TI, ', ex. path: ' EP, ', ex. distance: ' ED, ', rotation: ' RP]});
hold off; 

% save plot
file_name = ['Plot_' ID '_' Session '_' Trial '.jpeg'];
file_path = fullfile(folder, file_name);
saveas(wfig,file_path);

end
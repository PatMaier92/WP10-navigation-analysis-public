function [x_line, y_line, x_line_chosen, y_line_chosen, origin_x_line, origin_y_line,...
    x_line_ego, y_line_ego, goal_x_ego, goal_y_ego, ego_alley,...
    ideal_path, ideal_path_chosen, ideal_ego_path]=computeStartDependentIdealVariables(...
    Graph, graph_x, graph_y, start, goal, chosen_x, chosen_y,...
    alley_full_x, alley_full_y, rec_x, rec_y, cp_polyshape, polyshape_array)
% computeStartDependentIdealVariables: Function for determining start-dependent
% variables in Starmaze WP10. Requires Matlab 2021a for 'shortestpath' function and 
% 'interparc' function by John D'Errico (Matlab File Exchanger).
%
% Input:
% Graph is graph representation of all start-goal connections
% graph_x, graph_y are the x-/y-coordinates of the Graph
% start, goal are identifying integers
% chosen_x, chosen_y are the final x-/y-coordinates
% alley_full_*, rec_* are x-/y-coordinate vectors in polyshape-ready form 
% (i.e. repeating initial x/y for closed shape)
% cp_polyshape is inner ring as one combined polyshape
% polyshape_array is array of all polyshape elements
%
% Returns:
% *x_line*, *y_line* are vectors with ideal x-/y-coordinates
% goal_x_ego, goal_y_ego are x-/y-coordinates of the egocentric goal
% ego_alley is integer for hypothetical egocentric goal alley 
% ideal_path* are ideal path length values

%% shortest path from original start to goal
start_node=7; 
goal_node=size(Graph.Nodes,1)+1-goal;
[origin_x_line, origin_y_line, ~]=computeShortestPath(true,...
    cp_polyshape, polyshape_array, Graph, graph_x, graph_y, start_node, goal_node,...
    0, 0);

% % helper plot
% [path_nodes,~]=shortestpath(Graph, start_node, goal_node);
% figure; plot(polyshape_array); hold on; 
% pl = plot(Graph,'XData',graph_x,'YData',graph_y);
% highlight(pl,path_nodes,'EdgeColor','r');
% xlim([0 1]); ylim([0 1]); hold off; 
% clear path_nodes; 
% % test plot
% figure; plot(polyshape_array); hold on;
% plot(origin_x_line, origin_y_line, 'k-');
% xlim([0 1]); ylim([0 1]); hold off;

%% shortest path from actual start to egocentric goal (using a rotation matrix) 
% define x- and y-data for original line
v = [origin_x_line ; origin_y_line];
% define center of rotation
x_center = 0.5; y_center = 0.5;
% create a matrix
center = repmat([x_center; y_center], 1, length(v));

% define rotation matrix
if start==9
    theta=-360/5*1; % to rotate 72° clockwise
elseif start==1
     theta=-360/5*2; 
elseif start==3
     theta=-360/5*3; 
elseif start==5
     theta=-360/5*4; 
else
     theta=-360/5*0; % no rotation for original and inner starts
end
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

% do rotation
vo = R*(v - center) + center;

% get rotated x- and y-data
r_x_line = vo(1,:);
r_y_line = vo(2,:);

% correct rotated x- and y-data to account for measurement/rotation errors
% i.e., find nearest vertex in actual starmaze boundaries
% otherwise, your egocentric paths might be slightly off/outside the maze. 
[vid,~,~] = nearestvertex(cp_polyshape,r_x_line(2:end-1),r_y_line(2:end-1));

% save ego path and goal location
if ~mod(start,2) % dummy for inner starts (even start integer)
    x_line_ego=[999; 998]; y_line_ego=[999; 998]; 
    goal_x_ego=0; goal_y_ego=0; 
else 
    x_line_ego=[r_x_line(1); cp_polyshape.Vertices(vid,1); r_x_line(end)]; 
    y_line_ego=[r_y_line(1); cp_polyshape.Vertices(vid,2); r_y_line(end)]; 
    goal_x_ego=r_x_line(end); goal_y_ego=r_y_line(end);
end 

% % test plot
% figure; plot(cp_polyshape); hold on;
% plot(origin_x_line, origin_y_line, 'k-', x_line_ego, y_line_ego, 'rx', x_center, y_center, 'bo');
% xlim([0 1]); ylim([0 1]); hold off;

% get egocentric final alley integer
ego_alley=0;
[~,col]=size(alley_full_x);
for c=1:col
    if inpolygon(goal_x_ego,goal_y_ego,alley_full_x(:,c),alley_full_y(:,c)) 
        ego_alley=c*2-1;
    elseif inpolygon(goal_x_ego,goal_y_ego,rec_x(:,c),rec_y(:,c))
        ego_alley=c*2;
    end
end

% calculate ideal ego path length value (external function)
ideal_ego_path=computePathLength(x_line_ego, y_line_ego);

%% shortest path from actual start to goal
start_node=start;
goal_node=size(Graph.Nodes,1)+1-goal;
[x_line, y_line, ideal_path]=computeShortestPath(true,...
    cp_polyshape, polyshape_array, Graph, graph_x, graph_y, start_node, goal_node,...
    0, 0);

% % test plot
% [path_nodes,~]=shortestpath(Graph, start_node, end_node);
% figure; plot(polyshape_array); hold on; 
% pl = plot(Graph,'XData',graph_x,'YData',graph_y);
% highlight(pl,path_nodes,'EdgeColor','r');
% xlim([0 1]); ylim([0 1]); hold off; 
% clear path_nodes; 
% % test plot
% figure; plot(polyshape_array); hold on;
% plot(x_line, y_line, 'k-');
% xlim([0 1]); ylim([0 1]); hold off;

%% shortest path from actual start to chosen goal 
[x_line_chosen, y_line_chosen, ideal_path_chosen]=computeShortestPath(false,...
    cp_polyshape, polyshape_array, Graph, graph_x, graph_y, start, 0,...
    chosen_x, chosen_y);
    
end
function [x_line, y_line, ideal_length]=computeShortestPath(toGoal, ...
    central_polygon, full_polyshape, graph, graph_x, graph_y, ...
    start_node, goal_node, end_x, end_y)
% computeShortestPath: Helper function for computing the shortest path
% between any two locations in the WP10 Starmaze environment.
% Requires Matlab 2021a for 'shortestpath' function and 
% 'interparc' function by John D'Errico (Matlab File Exchanger).
% 
% Input:
% toGoal is boolean, indicating if path ends at a goal (i.e. ends at graph node)
% central_polygon is inner ring as one combined polyshape
% full_polyshape is full combined polyshape 
% graph is graph representation of all start-goal connections
% graph_x, graph_y are the x-/y-coordinates of the graph
% start_node is integer for start in graph
% goal_node is integer for goal in graph
% end_x, end_y are x-/y-coordinates of end point
% 
% Returns: 
% x_line, y_line are vectors with ideal x-/y-coordinates
% ideal_length ideal path length values

% boolean and integer indicating if end point is in graph
[inGraph,end_node]=ismembertol([end_x end_y],[graph_x' graph_y'],1e-9,'ByRows',true);

%% if end point of the path is part of the graph
if toGoal || inGraph
    % get nodes 
    if toGoal 
        [path_nodes,~]=shortestpath(graph, start_node, goal_node);
    elseif inGraph
        [path_nodes,~]=shortestpath(graph, start_node, end_node);
    end 

    % set path with all nodes
    x1=graph_x(path_nodes);
    y1=graph_y(path_nodes);
    length1=computePathLength(x1, y1);
    
    % set path with reduced nodes
    x2=graph_x([path_nodes(1:end-2) path_nodes(end)]);
    y2=graph_y([path_nodes(1:end-2) path_nodes(end)]);
    length2=computePathLength(x2, y2);
    % interpolate last segment (using 'interparc' function by John D'Errico (Matlab File Exchanger))
    [xi_s2,yi_s2]=interpolateData(x2(end-1:end),y2(end-1:end),length2);
    % check if last segment is valid (i.e., not outside starmaze area)
    isValid=all(isinterior(union(full_polyshape),xi_s2,yi_s2));
    
    % determine best option and set values
    if length2 < length1 && isValid
        x_line=x2; y_line=y2; ideal_length=length2; 
    else
        x_line=x1; y_line=y1; ideal_length=length1; 
    end
    
%% if end point of the path is NOT part of the graph
else   
    % find nearest vertices in central polygon
    % nearest
    [vid,~,~] = nearestvertex(central_polygon,end_x,end_y);
    node_x1=central_polygon.Vertices(vid,1);
    node_y1=central_polygon.Vertices(vid,2);
    % second nearest
    central_polygon.Vertices(vid,:)=[]; % remove nearest vertex and check again
    [vid,~,~] = nearestvertex(central_polygon,end_x,end_y);
    node_x2=central_polygon.Vertices(vid,1);
    node_y2=central_polygon.Vertices(vid,2); clear vid central_polygon;
    
    % get nodes
    % for nearest
    [~,node_index1]=ismembertol([node_x1 node_y1],[graph_x' graph_y'],1e-9,'ByRows',true);
    [path_nodes1,~]=shortestpath(graph,start_node,node_index1);
    % for second nearest
    [~,node_index2]=ismembertol([node_x2 node_y2],[graph_x' graph_y'],1e-9,'ByRows',true);
    [path_nodes2,~]=shortestpath(graph,start_node,node_index2); clear node*;
  
    % create four path options with added goal location x-/y-coordinates
    % 1a) all nodes to nearest vertex
    x1_all=[graph_x(path_nodes1) end_x];
    y1_all=[graph_y(path_nodes1) end_y];
    length1_all=computePathLength(x1_all, y1_all);
    
    % 1b) reduced nodes to nearest vertex
    x1_reduced=[graph_x(path_nodes1(1:end-1)) end_x];
    y1_reduced=[graph_y(path_nodes1(1:end-1)) end_y];
    length1_reduced=computePathLength(x1_reduced, y1_reduced);
    % interpolate last segment (using 'interparc' function by John D'Errico (Matlab File Exchanger))
    [xi_s1,yi_s1]=interpolateData(x1_reduced(end-1:end),y1_reduced(end-1:end),length1_reduced);
    % check if last segment is valid (i.e., not outside starmaze area)
    p1Valid=all(isinterior(union(full_polyshape),xi_s1,yi_s1)); clear *_s1;
    
    % 2a) all nodes to second nearest vertex
    x2_all=[graph_x(path_nodes2) end_x];
    y2_all=[graph_y(path_nodes2) end_y];
    length2_all=computePathLength(x2_all, y2_all);
    
    % 2b) reduced nodes to second nearest vertex
    x2_reduced=[graph_x(path_nodes2(1:end-1)) end_x];
    y2_reduced=[graph_y(path_nodes2(1:end-1)) end_y];
    length2_reduced=computePathLength(x2_reduced, y2_reduced);
    % interpolate last segment (using 'interparc' function by John D'Errico (Matlab File Exchanger))
    [xi_s2,yi_s2]=interpolateData(x2_reduced(end-1:end),y2_reduced(end-1:end),length2_reduced);
    % check if last segment is valid (i.e., not outside starmaze area)
    p2Valid=all(isinterior(union(full_polyshape),xi_s2,yi_s2)); clear *_s2;
    
    % determine best option and set values
    % compare 1a) and 1b)
    if length1_reduced < length1_all && p1Valid
        x_line1=x1_reduced; y_line1=y1_reduced; length1=length1_reduced;
    else
        x_line1=x1_all; y_line1=y1_all; length1=length1_all;
    end
    % compare 2a) and 2b)
    if length2_reduced < length2_all && p2Valid
        x_line2=x2_reduced; y_line2=y2_reduced; length2=length2_reduced;
    else
        x_line2=x2_all; y_line2=y2_all; length2=length2_all;
    end
    % compare 1) and 2)
    if length1 < length2
        x_line=x_line1; y_line=y_line1; ideal_length=length1;
    else
        x_line=x_line2; y_line=y_line2; ideal_length=length2;
    end
    clear length* *_p1_* *_p2_*;
end
 
% % % test plot
% % figure; plot(polyshape_array); hold on; 
% % plot(end_x, end_y, 'ro', x1_all, y1_all, 'm-', x2_all, y2_all, 'r-',...
% %     x1_reduced, y1_reduced,'b--', x2_reduced, y2_reduced, 'y--');
% % xlim([0 1]); ylim([0 1]); hold off; 

end
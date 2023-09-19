function plotTestFigure(session,polyshape,goal_x,goal_y,start_x,start_y,goal_locs,G,graph_x,graph_y)
% plotTestFigure Create Starmaze WP10 testfigure for visualization of all
% xy-coordinates, goals and starting points. 

% assign variables for practise or starmaze plot
if session=="s_maze"
    my_title="Starmaze WP10"; 
    start_index=7;
elseif session=="p_maze"
    my_title="Practise Maze WP10";
    start_index=1; 
end

% plot figure
figure('Position',[500 200 580 500]);
set(gca,'xtick',[0 1],'ytick',[0 1]);
plot(polyshape);
axis([0 1 0 1]);
title(my_title);
hold on
for g=1:length(goal_x)
    viscircles([goal_x(g) goal_y(g)], 0.01, 'Color','red');
end
viscircles([start_x(start_index) start_y(start_index)], 0.01, 'Color','blue');
if session=="s_maze"
    p=plot(G,'XData',graph_x,'YData',graph_y);
    labelnode(p,[1:28],G.Nodes.Names);
    text(0.1,0.75,'Alley 5 - I');
    text(0.8,0.75,'Alley 2 - C');
    text(0.05,0.2,'Alley 4 - G');
    text(0.8,0.2,'Alley 3 - E');
    text(0.6,0.9,'Alley 1 - A');
else 
    for i=1:length(goal_locs)
        text(goal_x(i)+0.02, goal_y(i), goal_locs(i))
    end 
end 
text(start_x(start_index)-0.05, start_y(start_index)-0.05, 'Initial start'); 
hold off

end
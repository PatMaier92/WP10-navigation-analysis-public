function [x_vec, y_vec]=computeRandomLocations(polyshape, n_points)
% computeRandomLocations Creates set of randomly and unformly distributed
% coordinates inside the Starmaze WP10. 
%
% Input: 
% polyshape is an array with the Starmaze polyshape 
% n_points is number of desired points
%
% Returns:  
% x_in, y_in are vectors with randomly and uniformly distributed x-/y-coordinates
    
% define parameters
m = 10000; % number of points
R = 0.5; % radius of a circle
x0 = 0.5; % center of a circle
y0 = 0.5; % center of a circle

% create random set of data points in a circle 
rng(10); t = 2*pi*rand(m,1); 
r = R*sqrt(rand(m,1));
x = x0 + r.*cos(t);
y = y0 + r.*sin(t);

% find those data points that overlap with the Starmaze shape
union_poly = union(polyshape);
i = inpolygon(x,y,union_poly.Vertices(:,1),union_poly.Vertices(:,2));
x_vec = x(i);
y_vec = y(i);

% reduce the random data points to n_points
n = numel(x(i));
j = sort(randperm(n,n_points));
x_vec = x_vec(j);
y_vec = y_vec(j);

% % helper plot
% figure; plot(polyshape); hold on;
% plot(x_vec,y_vec, 'ro', 'MarkerSize', 5);
% hold off;

end
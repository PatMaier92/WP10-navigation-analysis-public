function [tri_x,tri_y,tri]=setTriPolyshape(n_alleys,alley_x,alley_y,pentagon_x,pentagon_y)
% setTriPolyshape: Creates a polyshape of Starmaze WP10 triangles. 
%
% Input: n_alleys is integer with number of available alleys, alley_x and
% alley_y are arrays with x-/y-coordinates of alley corners, pentagon_y and
% pentagon_y are vectors with x-/y-coordinates of inner pentagon. 
%
% Returns: triangle polyshapes

for a=1:n_alleys
    tri_x(:,a)=[alley_x(4,a);alley_x(3,a);pentagon_x(1,a);alley_x(4,a)];
    tri_y(:,a)=[alley_y(4,a);alley_y(3,a);pentagon_y(1,a);alley_y(4,a)];
    tri{a}=polyshape(tri_x(:,a),tri_y(:,a));
end

end

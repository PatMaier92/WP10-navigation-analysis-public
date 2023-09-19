function [rec_x,rec_y,rec]=setRectPolyshape(n_alleys,alley_x,alley_y,pentagon_x,pentagon_y)
% setRectPolyshape: Creates a polyshape of Starmaze WP10 rectangles. 
% 
% Input: n_alleys is integer with number of available alleys, alley_x and
% alley_y are arrays with x-/y-coordinates of alley corners, pentagon_y and
% pentagon_y are vectors with x-/y-coordinates of inner pentagon. 
%
% Returns: rectangle polyshapes.

for a=1:n_alleys
    if a==n_alleys
        x=-n_alleys;
    else
        x=0;
    end
    
    rec_x(:,a)=[alley_x(3,a);alley_x(4,a+1+x);pentagon_x(1,a+1+x);pentagon_x(1,a);alley_x(3,a)];
    rec_y(:,a)=[alley_y(3,a);alley_y(4,a+1+x);pentagon_y(1,a+1+x);pentagon_y(1,a);alley_y(3,a)];
    rec{a}=polyshape(rec_x(:,a),rec_y(:,a));
end

end
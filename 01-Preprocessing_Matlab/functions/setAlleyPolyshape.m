function [alley_full_x,alley_full_y,alley_polyshape,...
    alley_half_out_x, alley_half_out_y, alley_polyshape_1,...
    alley_half_in_x, alley_half_in_y, alley_polyshape_2]=setAlleyPolyshape(alley_x,alley_y)
% setAlleyPolyshape: Creates a polyshape of Starmaze WP10 alleys. 
%
% Input: alley_x, alley_y are 4*5 arrays with x and y coordinates corners.
%
% Returns: Polyshapes for all alleys, both full and half alleys and combined. 

[~, n]=size(alley_x);

for alley=1:n
    alley_full_x(:,alley)=[alley_x(1,alley);alley_x(2,alley);alley_x(3,alley);alley_x(4,alley);alley_x(1,alley)];
    alley_full_y(:,alley)=[alley_y(1,alley);alley_y(2,alley);alley_y(3,alley);alley_y(4,alley);alley_y(1,alley)];
    alley_polyshape{alley}=polyshape(alley_full_x(:,alley),alley_full_y(:,alley));
    
    x1=alley_x(1,alley)+0.5*(alley_x(4,alley)-alley_x(1,alley));
    x2=alley_x(2,alley)+0.5*(alley_x(3,alley)-alley_x(2,alley));
    y1=alley_y(1,alley)+0.5*(alley_y(4,alley)-alley_y(1,alley));
    y2=alley_y(2,alley)+0.5*(alley_y(3,alley)-alley_y(2,alley));
    
    alley_half_out_x(:,alley)=[alley_x(1,alley); alley_x(2,alley); x2;x1; alley_x(1,alley)];
    alley_half_out_y(:,alley)=[alley_y(1,alley); alley_y(2,alley); y2;y1; alley_y(1,alley)];
    alley_polyshape_1{alley}=polyshape(alley_half_out_x(:,alley),alley_half_out_y(:,alley));
    
    alley_half_in_x(:,alley)=[alley_x(3,alley); alley_x(4,alley); x1;x2; alley_x(3,alley)];
    alley_half_in_y(:,alley)=[alley_y(3,alley); alley_y(4,alley); y1;y2; alley_y(3,alley)];
    alley_polyshape_2{:,alley}=polyshape(alley_half_in_x(:,alley),alley_half_in_y(:,alley));
end
end
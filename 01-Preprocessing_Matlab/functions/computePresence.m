function [abs_presence, rel_presence, abs_time_in_zone, index]=computePresence(alley_int,...
    x, y, alley_poly, optional_tri_poly, time)
% computePresence: Compute absolute and relative presence in zone. 
% 
% Input: 
% alley_int is zone identifier (integer) 
% x, y are x-/y-coordinates (float) 
% alley_poly, tri_poly are polyshapes
% time is total time (float)
% 
% Returns:
% abs_presence (integer), rel_presence (float), 
% abs_time_in_zone(float), index (boolean vector)

% convert alley_int (from [1 3 5 7 9] to [1 2 3 4 5])
alley_int=(alley_int+1)/2;

% get index(es) in zone 
target_poly=union(alley_poly{alley_int}, optional_tri_poly{alley_int});
index=inpolygon(x,y,target_poly.Vertices(:,1),target_poly.Vertices(:,2));
% figure; plot(sm.coord.full_poly); hold on;
% plot(target_poly, 'FaceColor', 'k');
% plot(x(index), y(index), 'rx');

% compute values
abs_presence=numel(x(index));
rel_presence=abs_presence/length(x);
abs_time_in_zone=time*rel_presence;
                
end

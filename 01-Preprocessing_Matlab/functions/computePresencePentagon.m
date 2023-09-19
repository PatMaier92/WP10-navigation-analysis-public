function [abs_presence, rel_presence, abs_time_in_zone, index]=computePresencePentagon(alley_int,...
    x, y, rec_poly, tri_poly, time)
% computePresencePentagon: Compute absolute and relative presence in zone.
% This is a variant of computePresence, specifically for inner zones. 
% 
% Input: 
% alley_int is zone identifier (integer) 
% x, y are x-/y-coordinates (float) 
% rec_poly, tri_poly are polyshapes
% time is total time (float)
% 
% Returns:
% abs_presence (integer), rel_presence (float), 
% abs_time_in_zone(float), index (boolean vector)

% convert alley_int (from [2 4 6 8 10] to [1 2 3 4 5])
alley_int=(alley_int)/2;

% get index(es) in zone 
if alley_int < 5 
    target_poly=union(tri_poly{alley_int}, tri_poly{alley_int+1});
elseif alley_int == 5
    target_poly=union(tri_poly{5}, tri_poly{1});
end 
target_poly=union(target_poly, rec_poly{alley_int});  
index=inpolygon(x,y,target_poly.Vertices(:,1),target_poly.Vertices(:,2));
% figure; plot(sm.coord.full_poly); hold on;
% plot(target_poly, 'FaceColor', 'k');
% plot(x(index), y(index), 'rx');

% compute values
abs_presence=numel(x(index));
rel_presence=abs_presence/length(x);
abs_time_in_zone=time*rel_presence;
                
end

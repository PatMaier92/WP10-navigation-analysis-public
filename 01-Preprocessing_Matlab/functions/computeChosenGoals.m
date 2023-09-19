function [goal_i, alley_s, alley_i, obj_at_location]=computeChosenGoals(rand_dict, goal_s, ...
    alley_poly, rec_poly, tri_poly, alley_names, goal_names, xend, yend)
% computeChosenGoals: Returns chosen goal location as integer and string.  
%
% Input: 
% rand_dict with randomization info (structure/dictionary)
% goal_s (string) only relevant for "Timeout". The goal string 
%   is error-prone for triangles/intersections and therefore recalculated by this function. 
% alley_poly, rec_poly, tri_poly contain polyshapes
% alley_names, goal_names are string vectors
% xend, yend are chosen final x-/y-coordinates
%
% Returns: 
% goal_i (integer), alley_s (string), alley_i (integer), obj_at_location (string) 

if goal_s == "Timeout"
    alley_i = 999; 
    alley_s = "Timeout"; 
    goal_i = 999; 
    obj_at_location = "999";  
    disp('Timeout in computeChosenGoals.m. Set to 999.\n');
else 
    for c=1:5
        if inpolygon(xend, yend, alley_poly{c}.Vertices(:,1), alley_poly{c}.Vertices(:,2)) ...
                || inpolygon(xend, yend, tri_poly{c}.Vertices(:,1), tri_poly{c}.Vertices(:,2)) % chosen location in alley or triangle (outer arm)
            % chosen alley integer
            alley_i = c*2-1;
            break; 
        elseif inpolygon(xend, yend, rec_poly{c}.Vertices(:,1), rec_poly{c}.Vertices(:,2)) % chosen location in rectangle (inner arm)
             % chosen alley integer
            alley_i = c*2;
            break;
        end 
    end 
    
    % chosen alley string
    alley_s = alley_names(alley_i);
    
    % chosen goal string
    goal_i = find(contains(goal_names, alley_s)); 
    if isempty(goal_i)
        goal_i = 999;
        %fprintf('Unknown goal int input %s in computeChosenGoals.m. Set to 999.\n', goal_s);
    end 

    % correct object at chosen location 
    obj_at_location="999"; 
    fields = fieldnames(rand_dict); 
    for i=1:length(fields)
        if strcat("M", alley_s) == string(fields{i})
            key = fields{i};
            obj_at_location = string(rand_dict.(key).object);
        end
    end 
end 

end

function memory_score=computeMemoryScore(distribution, final_distance,...
    goal_int, alley_int)
% computeMemoryScore: Computes memory score, i.e. final distance to goal 
% is set in relation to randomly distributed points in the WP10 Starmaze. 
%
% Input: 
% distribution is matrix of final distance values to random points (float)
% final_distance is final distance value for this trial (float)
% goal_int indicates goal location for this comparison (integer) 
% alley_int indicates desired alley for this comparison (integer)
%
% Returns: standardized memory_score (0=low, 1=high) (float) 

% convert alley_int (from [1 3 5 7 9] to [1 2 3 4 5])
alley_int=(alley_int+1)/2; 

% compute score 
percentage_below = sum(distribution{goal_int,alley_int} < final_distance) / 1000;
memory_score = 1 - percentage_below;
  
end


function [sum_rotation]=computeRotation(r)
% computeRotation: Calculates total rotation in degrees as sum of
% absolute change in yaw rotation (first derivative). This value includes 
% rotation due to x-/y-trajectory (i.e. left-forward movement).
%
% Input: 
% r is vector with unwrapped (!) z-coordinates with yaw rotation (float) 
% 
% Returns:
% sum_rotation (float)

sum_rotation=sum(abs(diff(r)));
                
end

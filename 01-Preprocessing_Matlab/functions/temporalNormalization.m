function [ti, xi, yi, ri]=temporalNormalization(nsr, t, x, y, r)
% temporalNormalization: Interpolation of coordinates. 

% Input: 
% nsr is desired new sampling rate.
% t, x, y, r are original data vectors.
%
% Returns: ti, xi, yi, ri are interpolated data vectors.  

% compute equally spaced time vector (query points) 
ti=[t(1):nsr:t(end)]';

% compute spatio-temporal interpolation
xi=[interp1(t,x,ti,'linear'); x(end)];
yi=[interp1(t,y,ti,'linear'); y(end)];
ri=[interp1(t,r,ti,'linear'); r(end)];

% add end point
ti=[ti; t(end)]; 

% % test plot xy
% figure; plot(sm.coord.full_poly); hold on;
% plot(x, y, 'k+');
% plot(xi, yi, 'r*');
% xlim([0 1]); ylim([0 1]); hold off;
% 
% % test plot r
% figure; hold on;
% plot(r, 'k+');
% plot(ri, 'r+');
% hold off;

end 



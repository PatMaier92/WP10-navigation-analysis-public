function DN=setNormalizedValues(c,cmin,cmax)
% setNormalizedValues: Normalizes coordinates based on min/max values. 

% Input: 
% c is input value to be normalized.
% cmin, cmax are minumum, maximum values.
%
% Returns: DN is normalized value. 

DN=((c-cmin)/(cmax-cmin));



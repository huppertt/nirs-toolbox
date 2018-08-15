function prop = brain(lambda, so2, hbt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 2, so2 = 0.7; end
    if nargin < 3, hbt = 60; end
    if nargin < 1, lambda = [690 830]; end
    
    prop = nirs.media.SpectralProp( so2, hbt, lambda );

end


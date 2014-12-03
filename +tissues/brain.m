function prop = brain( so2, hbt, lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 1, so2 = 0.7; end
    if nargin < 2, hbt = 60; end
    if nargin < 3, lambda = [690 830]; end
    
    prop = nirs2.SpectralProperties( so2, hbt, lambda );

end


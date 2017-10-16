function [R,p,dfe]=grangers_compat(varargin)
% [R,p,dfe]=grangers_compat(data,modelorder,robust_flag,includeZeroLag)
%
% This is a wrapper for the grangers code to maintain compatibility with
%   the connectivity modules and sFC data class.
% The F-stats are converted to R-values, which will likely need to be 
%   centered by subtracting a null condition

[G,F,dfe1,dfe2,p] = nirs.sFC.grangers(varargin{:});

% Remove directionality
F = max( F , permute(F,[2 1 3]) );

% Convert to Fisher Z-distribution
Z = log( F ) / 2;

% Convert to r-value
R = tanh(Z);

assert(isreal(R),'r-values should not be complex!');

dfe = nanmean(dfe1(:)) + 1i * nanmean(dfe2(:));

end
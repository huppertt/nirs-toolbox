%NORMALIZE Normalize coefficients between -1 and 1.
%   G = NORMALIZE(Hd) normalizes the feed-forward coefficients between -1
%   and 1 and returns the gain G due to normalization.  Subsequent calls to
%   NORMALIZE will not change the feed-forward coefficients and G will
%   always return the gain used in the first normalization.
%
%   See also DFILT/DENORMALIZE.
% Copyright 1988-2004 The MathWorks, Inc.

% [EOF]

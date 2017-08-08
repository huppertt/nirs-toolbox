function y = downsample(x,N,varargin)
%DOWNSAMPLE Downsample input signal.
%   DOWNSAMPLE(X,N) downsamples input signal X by keeping every
%   N-th sample starting with the first. If X is a matrix, the
%   downsampling is done along the columns of X.
%
%   DOWNSAMPLE(X,N,PHASE) specifies an optional sample offset.
%   PHASE must be an integer in the range [0, N-1].
%
%   % Example 1:
%   %   Decrease the sampling rate of a sequence by 3.
%
%   x = [1 2 3 4 5 6 7 8 9 10];
%   y = downsample(x,3)
%
%   % Example 2:
%   %   Decrease the sampling rate of the sequence by 3 and add a 
%   %   phase offset of 2.
%
%   x = [1 2 3 4 5 6 7 8 9 10];
%   y = downsample(x,3,2)
%
%   % Example 3:
%   %   Decrease the sampling rate of a matrix by 3.
%
%   x = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
%   y = downsample(x,3)
%
%   See also UPSAMPLE, UPFIRDN, INTERP, DECIMATE, RESAMPLE.

%   Copyright 1988-2002 The MathWorks, Inc.

y = updownsample(x,N,'Down',varargin{:});

% [EOF] 

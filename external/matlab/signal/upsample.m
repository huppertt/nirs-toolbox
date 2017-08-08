function y = upsample(x,N,varargin)
%UPSAMPLE Upsample input signal.
%   UPSAMPLE(X,N) upsamples input signal X by inserting
%   N-1 zeros between input samples.  X may be a vector
%   or a signal matrix (one signal per column).
%
%   UPSAMPLE(X,N,PHASE) specifies an optional sample offset.
%   PHASE must be an integer in the range [0, N-1].
%
%   % Example 1:
%   %   Increase the sampling rate of a sequence by 3.
%
%   x = [1 2 3 4];      % Defining data
%   y = upsample(x,3)   % Upsample input signal            
%
%   % Example 2:
%   %   Increase the sampling rate of the sequence by 3 and add a   
%   %   phase offset of 2.
%
%   x = [1 2 3 4];      % Defining data
%   y = upsample(x,3,2) % Upsammple by 3 and adding phase offset of 2 
%
%   % Example 3:
%   %   Increase the sampling rate of a matrix by 3. 
%
%   x = [1 2; 3 4; 5 6;];   % Defining data
%   y = upsample(x,3)       % Increasing sampling rate      
%
%   See also DOWNSAMPLE, UPFIRDN, INTERP, DECIMATE, RESAMPLE.

%   Copyright 1988-2010 The MathWorks, Inc.

y = updownsample(x,N,'Up',varargin{:});

% [EOF] 

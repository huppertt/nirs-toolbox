function w = blackman(varargin)
%BLACKMAN   Blackman window.
%   BLACKMAN(N) returns the N-point symmetric Blackman window in a column
%   vector.
%   BLACKMAN(N,SFLAG) generates the N-point Blackman window using SFLAG
%   window sampling. SFLAG may be either 'symmetric' or 'periodic'. By 
%   default, a symmetric window is returned. 
%
%   % Example:
%   %   Create a 64-point Blackman window and display the result using 
%   %   WVTool.
%
%   L=64;
%   wvtool(blackman(L))
%
%   See also  HAMMING, HANN, WINDOW.

%   Copyright 1988-2013 The MathWorks, Inc.

% Check number of inputs
narginchk(1,2);

[w,msg,msgobj] = gencoswin('blackman',varargin{:});
if ~isempty(msg), error(msgobj); end

% [EOF] blackman.m

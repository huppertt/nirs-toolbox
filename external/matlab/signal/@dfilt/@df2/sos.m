function Hsos = sos(Hd,varargin)
%SOS  Convert to second-order-sections.
%   Hsos = SOS(Hd) converts discrete-time filter Hd to second-order section
%   form.
%
%   For optional arguments SOS(Hd,...), see TF2SOS.
%
%   Example:
%          [b,a] = butter(8,.5);
%          Hd = dfilt.df2(b,a);
%          Hsos = sos(Hd,'up',inf)
%  
%   See also DFILT.

%   Author: Thomas A. Bryan, V. Pellissier
%   Copyright 1988-2008 The MathWorks, Inc.

% Convert to second-order-section matrix and gain.
[b,a] = tf(Hd);
[s,g] = tf2sos(b,a,varargin{:});

% Build a cell array of all the numerators and denominators from the sos
% array. 
n = size(s,1);
for k=1:n,
  c{2*k-1} = s(k,1:3);
  c{2*k}   = s(k,4:6);
end
c{end+1} = g;

% Cascade in all the numerators and denominators using the same filter
% structure as this Hd.
Hsos = thissos(Hd,c);

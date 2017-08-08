function Hsos = sos(Hd,varargin)
%SOS  Convert to second-order-sections.
%   Hsos = SOS(Hd) converts discrete-time filter Hd to second-order section
%   form.
%
%   SOS(Hd,DIR_FLAG) specifies the ordering of the 2nd order sections. If
%   DIR_FLAG is equal to 'UP', the first row will contain the poles closest
%   to the origin, and the last row will contain the poles closest to the
%   unit circle. If DIR_FLAG is equal to 'DOWN', the sections are ordered
%   in the opposite direction. The zeros are always paired with the poles
%   closest to them. DIR_FLAG defaults to 'UP'.
%
%   EXAMPLE:
%          [b,a] = butter(8,.5);
%          Hd = dfilt.df2(b,a);
%          Hsos = sos(Hd,'up',inf)
%  
%   See also DFILT.

%   Author: Thomas A. Bryan
%   Copyright 1988-2008 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct')); % No scaling allowed

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

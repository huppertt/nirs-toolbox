function y = sawtooth(t,width)
%SAWTOOTH Sawtooth and triangle wave generation.
%   SAWTOOTH(T) generates a sawtooth wave with period 2*pi for the
%   elements of time vector T.  SAWTOOTH(T) is like SIN(T), only
%   it creates a sawtooth wave with peaks of +1 to -1 instead of
%   a sine wave.
%
%   SAWTOOTH(T,WIDTH) generates a modified triangle wave where WIDTH, a
%   scalar parameter between 0 and 1, determines the fraction between 0
%   and 2*pi at which the maximum occurs. The function increases from -1
%   to 1 on the interval 0 to WIDTH*2*pi, then decreases linearly from 1
%   back to -1 on the interval WIDTH*2*pi to 2*pi. Thus WIDTH = .5 gives
%   you a triangle wave, symmetric about time instant pi with peak amplitude
%   of one.  SAWTOOTH(T,1) is equivalent to SAWTOOTH(T).
%
%   Caution: this function is inaccurate for huge numerical inputs
%
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   SINC and TRIPULS.

%   Copyright 1988-2012 The MathWorks, Inc.

if nargin == 1,
    width = 1;
end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(t,width);
catch ME
    throwAsCaller(ME);
end
  
if (width > 1) || (width < 0),
    error(message('signal:sawtooth:InvalidRange'))
end

rt = rem(t,2*pi)*(1/2/pi);
i1 = find( ((rt<=width)&(rt>=0)) | ((rt<width-1)&(rt<0)) );
i2 = 1:length(t(:));
i2(i1) = [];      % complement set
y = zeros(size(t));
y(i1) = ( ((t(i1)<0)&(rt(i1)~=0)) + rt(i1) - .5*width)*2;
if (width ~= 0),
    y(i1) = y(i1)*(1/width);
end
y(i2) = ( -(t(i2)<0) - rt(i2) + 1 - .5*(1-width))*2;
if (width ~= 1),
    y(i2) = y(i2)*(1/(1-width));
end
if width == 0
  y(rt == 0) = 1;
end


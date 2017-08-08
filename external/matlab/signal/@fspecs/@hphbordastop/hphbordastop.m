function h = hphbordastop
%HPHBORDASTOP   Construct a HPHBORDASTOP object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.hphbordastop;

h.ResponseType = 'Highpass halfband with filter order and stopband attenuation';
% [EOF]

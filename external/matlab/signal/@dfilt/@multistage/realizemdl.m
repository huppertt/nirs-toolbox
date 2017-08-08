function realizemdl(H,varargin)
%REALIZEMDL Filter realization (Simulink diagram).
%     REALIZEMDL(Hd) automatically generates architecture model of filter
%     Hd in a Simulink subsystem block using individual sum, gain, and
%     delay blocks, according to user-defined specifications.
%
%     REALIZEMDL(Hd, PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...) generates
%     the model with parameter/value pairs.
%
%    EXAMPLES:
%    % Realize multistage filter with map the coefficients to the ports
%    Hd = dfilt.dffir([0.05 0.9 0.05]);
%    Hgain = dfilt.scalar(2);
%    Hcas = dfilt.cascade(Hgain,Hd)
%    coeffnames.Stage1 = {'MyGain'};      % Use structure array to store the name
%    coeffnames.Stage2 = {'MyNum'};
%    realizemdl(Hcas,'MapCoeffsToPorts','on','CoeffNames',coeffnames);
%
%    See also DFILT/REALIZEMDL

%   Copyright 2004-2009 The MathWorks, Inc.

super_realizemdl_composite(H,varargin{:});

% [EOF]

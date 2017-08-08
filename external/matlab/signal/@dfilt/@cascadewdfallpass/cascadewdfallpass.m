function this = cascadewdfallpass(varargin)
%CASCADEWDFALLPASS Cascade of wave digital allpass structures.
%   Hd = DFILT.CASCADEWDFALLPASS(C1,C2,...) constructs a cascade of wave
%   digital allpass structures given the allpass coefficients in vectors
%   C1, C2, etc.
%
%   Each vector Ci represents one section of the cascade and must be a
%   vector of length 1, 2, or 4. When the length of each vector is Ci is 4,
%   the first and third coefficients in that vector must be zero. Each
%   section is a DFILT.WDFALLPASS filter with coefficients given by the
%   corresponding Ci vector.
%
%   Note that the vectors Ci do not have to be of the same length. For
%   example it is possible to cascade several fourth-order sections with
%   second-order, or first-order sections.
%
%   See the help for DFILT.WDFALLPASS for more information about the
%   vectors Ci and the transfer function of each section.
%
%   Note that one usually does not construct these filters directly.
%   Instead, they are used as part of the result of a design of an IIR
%   filter. See example below for more on this.
%
%   % Example #1: Design an IIR halfband filter that uses lattice wave
%   % digital filters, each branch of the parallel connection is a cascade
%   % wave digital allpass filter.
%   TW = 100;  % Transition width of filter to be designed, 100 Hz
%   Ast = 80;  % Stopband attenuation of filter to be designed, 80 db
%   Fs = 2000; % Sampling frequency of signal to be filtered
%   f = fdesign.halfband('TW,Ast',TW,Ast,Fs); % Store halfband design specs
%   % Now perform actual design. Hd contains two dfilt.cascadewdfallpass
%   Hd = design(f,'ellip','FilterStructure','cascadewdfallpass'); 
%   Hd.Stage(1).Stage(1) % Summary on one of the dfilt.cascadewdfallpass
%   realizemdl(Hd.stage(1).Stage(1)) % Requires Simulink; build model for filter
%
%   % Example #2: Construct a DFILT.CASCADEWDFALLPASS directly given
%   % allpass coefficients.
%   section1 = 0.8;
%   section2 = [1.5,0.7];
%   section3 = [1.8,0.9];
%   Hd = dfilt.cascadewdfallpass(section1,section2,section3);
%   info(Hd)    % Show information on the filter
%   fvtool(Hd)  % Visualize filter 
%
%   See also DFILT/STRUCTURES

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

this = dfilt.cascadewdfallpass;

this.FilterStructure = 'Cascade Wave Digital Filter Allpass';

constr(this,varargin{:});

% [EOF]

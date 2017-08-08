function this = fracdelay(varargin)
%FRACDELAY   Construct a fractional delay filter designer.
%   D = FDESIGN.FRACDELAY(DELTA) constructs a fractional delay filter
%   designer D with a delay of DELTA samples. The fractional delay DELTA
%   must be between 0 and 1.
%
%   D = FDESIGN.FRACDELAY(DELTA,'N') initializes the filter designer
%   'Specification' property to 'N' where the filter order N defaults to 3.
%
%   D = FDESIGN.FRACDELAY(DELTA,'N',N) initializes the filter designer
%   specifications with 'N' and sets the filter order to the value N. 
%
%   D = FDESIGN.FRACDELAY(DELTA,N) uses the default specification ('N') and
%   sets the filter order to the value N.
%
%   D = FDESIGN.FRACDELAY(...,Fs) specifies the sampling frequency (in Hz).
%   In this case, the fractional delay DELTA is expressed in seconds. The
%   fractional delay DELTA must be between 0 and 1/Fs.
%
%   % Example #1 
%        %Design a second order fractional delay filter of 0.2 samples using 
%        %the Lagrange method. Implement the filter using a Farrow FD structure.
%        d = fdesign.fracdelay(0.2,'N',2);
%        Hd = design(d, 'lagrange', 'FilterStructure', 'farrowfd'); 
%        fvtool(Hd, 'Analysis', 'grpdelay')
%
%   % Example #2
%        %Design a cubic fractional delay filter with a sampling frequency 
%        %of 8 kHz and a fractional delay of 50 microseconds using the 
%        %Lagrange method.
%        d = fdesign.fracdelay(50e-6,'N',3,8000);
%        Hd = design(d, 'lagrange', 'FilterStructure', 'farrowfd');
%        fvtool(Hd)
%
%   For more information about fractional delay filter implementations, see
%   the <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\farrowdemo.html'])">Fractional Delay Filters Using Farrow Structures</a> demo. 
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Author(s): V. Pellissier
%   Copyright 2005-2010 The MathWorks, Inc.

this = fdesign.fracdelay;
set(this, 'Response', 'Fractional Delay');
fd = 0.5;
if nargin>0,
    fd = varargin{1};
    varargin(1) = [];
end
setspecs(this, varargin{:});
this.FracDelay = fd;

% [EOF]

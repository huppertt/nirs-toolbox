function this = farrowfd(varargin)
%FARROWFD Fractional Delay Farrow filter.
%   Hd = DFILT.FARROWFD(D, COEFFS) constructs a discrete-time fractional
%   delay Farrow filter with COEFFS coefficients and D delay.
%
%   Farrow filters can be designed with the <a href="matlab:help fdesign.fracdelay">fdesign.fracdelay</a> filter designer. 
%
%   % EXAMPLE #1
%   coeffs = [-1/6 1/2 -1/3 0;1/2 -1 -1/2 1; -1/2 1/2 1 0;1/6 0 -1/6 0];
%   Hd = dfilt.farrowfd(.5, coeffs)
%   y = filter(Hd,1:10)
%
%   % EXAMPLE #2: Design a cubic fractional delay filter with the Lagrange method
%   fdelay = .2; % Fractional delay
%   d = fdesign.fracdelay(fdelay,'N',3);
%   Hd = design(d, 'lagrange', 'FilterStructure', 'farrowfd');
%   fvtool(Hd, 'Analysis', 'grpdelay') 
%
%   For more information about fractional delay filter implementations, see
%   the <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\farrowdemo.html'])">Fractional Delay Filters Using Farrow Structures</a> demo. 
%
%   See also DFILT

%   Copyright 2007-2010 The MathWorks, Inc.

error(nargchk(0,2,nargin,'struct'));

this = dfilt.farrowfd;
this.FilterStructure = 'Farrow Fractional Delay';
if nargin>0,
    this.FracDelay = varargin{1};
end
if nargin>1,
    this.Coefficients =  varargin{2};
end

% [EOF]

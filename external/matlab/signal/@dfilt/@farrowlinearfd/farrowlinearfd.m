function this = farrowlinearfd(varargin)
%FARROWLINEARFD Farrow Linear Fractional Delay filter.
%   Hd = DFILT.FARROWLINEARFD(D) constructs a discrete-time linear
%   fractional delay Farrow filter with the delay D.
%
%   % EXAMPLE
%   Hd = dfilt.farrowlinearfd(.5)
%   y = filter(Hd,1:10)
%
%   For more information about fractional delay filter implementations, see
%   the <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\farrowdemo.html'])">Fractional Delay Filters Using Farrow Structures</a> demo. 
%
%   See also DFILT

%   Copyright 2007-2010 The MathWorks, Inc.

error(nargchk(0,1,nargin,'struct'));
this = dfilt.farrowlinearfd;
this.FilterStructure = 'Farrow Linear Fractional Delay';
this.Coefficients = [-1 1;1 0];
this.States = 0;
if nargin==1,
    this.FracDelay = varargin{1};
end

% [EOF]

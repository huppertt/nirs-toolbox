function varargout = drawmask(this, hfm, hax, varargin)
%DRAWMASK   Draw the mask.

%   Copyright 2005-2012 The MathWorks, Inc.

narginchk(2,7);

% Make sure that we have a valid fmethod object.
hfm = getfmethod(this, hfm);

% Don't require a passed axes handle.  This is for demo purposes and should
% not be used by callers.
if nargin < 3
    hax = gca;
end

% Set up the mask utilities.
fcns = maskutils(this, isconstrained(hfm), hax, varargin{:});

% Get the mask information from the subclass.
[f, a] = getmask(this, fcns, 1); %Mask information assuming RCF of 1.
units = feval(fcns.getunits);
a = postprocessmask(hfm, a, units);

% avoid zeros on mask frequency vector since they will not be painted in
% log-scale
freqScale = feval(fcns.getfreqscale);
if strcmp(freqScale,'Log')
  xlim = feval(fcns.getxlim);
  f(f==0) = xlim(1); 
end

if nargout == 2
    varargout = {f, a};
else
    
    if isempty(hax)
        h = [];
    else
        % Put up the mask line.
        h = line(f, a, ...
            'Parent',    hax, ...
            'Color',     [0.8471 0.1608 0], ...
            'LineStyle', '--', ...
            'Tag',       'specification_mask');
    end
    
    if nargout == 1
        varargout = {h};
    else
        varargout = {f, a, h, fcns};
    end
end

% [EOF]

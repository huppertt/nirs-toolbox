function [xlim, ylim] = passbandzoom(this, hfm, varargin)
%PASSBANDZOOM   Returns the limits of the passband zoom.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% The method name is required, but everything else is optional.
error(nargchk(2,6,nargin,'struct'));

% If the method was passed in as a name, construct the appropriate object.
if ischar(hfm)
    hfm = feval(getdesignobj(this.CurrentSpecs, hfm));
end

% If we are passed a DFilt, use it.  Otherwise pass [] to the subclasses.
if nargin > 2 && isa(varargin{1}, 'dfilt.basefilter')
    Hd = varargin{1};
    varargin(1) = [];
else
    Hd = [];
end

% Call MASKUTILS to set up the nested function with the inputs which are
% the same as the DRAWMASK method.
fcns = maskutils(this, isconstrained(hfm), varargin{:});

% Call the subclasses THISPASSBANDZOOM to format the specs/measurements.
[xlim, ylim] = thispassbandzoom(this, fcns, Hd, hfm);

% The "ylimit" padding can be abstracted up but the "xlimit" cannot.
dr_ylim = diff(ylim);

% Add some buffer to the ylimits.
ylim = [ylim(1)-dr_ylim/10 ylim(2)+dr_ylim/10];

% [EOF]

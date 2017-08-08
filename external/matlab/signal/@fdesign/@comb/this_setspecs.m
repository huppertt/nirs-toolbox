function this_setspecs(this, varargin)
%THIS_SETSPECS   Set all the specs.
%   THIS_SETSPECS(S1, S2) Set the specifications in the order that they
%   appear in the 'SpecificationType' property.
%
%   THIS_SETSPECS(S1, S2, Fs) Set the Fs.  You must specify the Fs after all of
%   the other specifications for the current SpecificationType have been
%   specified, see example 2.
%
%   % Example #1:
%   h = fdesign.lowpass('n,fc')
%   h.setspecs(20, .4);
%   h
%
%   % Example #2:
%   h = fdesign.lowpass('n,fc')
%   h.setspecs(20, 4, 20);
%   h

%   Copyright 2008 The MathWorks, Inc.

if nargin > 2 && ischar(varargin{1}),
    set(this, 'SpecificationType', varargin{1});
    varargin(1) = [];
end

setspecs(this.CurrentSpecs, varargin{:});

% [EOF]

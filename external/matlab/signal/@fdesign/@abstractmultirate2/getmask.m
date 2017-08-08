function [F, A] = getmask(this, fcns, rcf, varargin)
%GETMASK   Get the mask.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

rcf  = getratechangefactors(this);

% Get the mask from the contained FDESIGN.
[F, A] = getmask(this.CurrentFDesign, fcns, max(rcf),varargin{:});

% [EOF]

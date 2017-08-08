function varargout = cumsec(this, indices, secondary)
%CUMSEC   Returns a vector of filters for the cumulative sections.
%   H = CUMSEC(Hd) returns a vector of SOS filter objects with the
%   cumulative sections.
%
%   H = CUMSEC(Hd, INDICES) return a vector of SOS filter objects whose
%   indices into the original filter are in INDICES.
%
%   H = CUMSEC(Hd, INDICES, SECONDARY) uses the secondary scaling points to
%   determine where the sections should be split when SECONDARY is true.
%   SECONDARY is false by default.  This option only has an effect on
%   DF2SOS and DF1TSOS filter objects. For these structures, the secondary
%   scaling points refers to the location between the recursive and the
%   nonrecursive part (i.e. the "middle" of the section).
%
%   CUMSEC(Hd,...) with no output arguments plots the magnitude response of
%   the cumulative sections using FVTOOL.
%
%   See also DFILT/SCALE, DFILT/SCALECHECK.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2,
    indices = 1:nsections(this);
    secondary = false;
elseif length(indices) == 1 & islogical(indices)
    secondary = indices;
    indices   = 1:nsections(this);
elseif nargin < 3,
    secondary = false;
end

if any(indices > nsections(this)),
    error(message('signal:dfilt:abstractsos:cumsec:exceedssecs'));
end

if nargout == 0,
    hopts = dspopts.sosview('View','Cumulative', 'SecondaryScaling', secondary);
    fvtool(this, 'SOSView', hopts);
else
    
    % If Secondary scaling points is requested and the structure supports them,
    % then move the numerator up.
    if secondary & shiftsecondary(this),

        % Make a copy of the filter so that we do not disturb the original
        % filters SOSMatrix.
        this = copy(this);

        sosM = this.sosMatrix;

        sosM(2:end, 1:3) = sosM(1:end-1, 1:3);
        sosM(1, 1:3) = [1 0 0];

        this.sosMatrix = sosM;
    end

    % Create an indexing variable because we don't always loop from 1:n.
    indx = 1;

    % Build copies and remove the extra sections with reorder.
    for jndx = indices
        h(indx) = copy(this);
        reorder(h(indx), 1:jndx);
        h(indx).ScaleValues(end) = 1;
        indx = indx+1;
    end
    varargout = {h};
end

% [EOF]

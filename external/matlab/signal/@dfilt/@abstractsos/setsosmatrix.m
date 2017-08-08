function s = setsosmatrix(this,s)
%SETSOSMATRIX Set the SOS matrix.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% Set the reference and check datatype
set(this, 'refsosMatrix', s);

% Set the new number of sections
oldnsections = this.nsections;
nsections = size(s,1);
this.nsections = nsections;

set_ncoeffs(this.filterquantizer, 6*nsections);

% Quantize the coefficients
quantizecoeffs(this);

if nsections~=oldnsections,
    % Reset the filter
    reset(this);
    set(this, 'ScaleValues', this.ScaleValues(1:end-this.NumAddedSV));
end

% Don't duplicate storage
s = []; 

clearmetadata(this);

% [EOF]

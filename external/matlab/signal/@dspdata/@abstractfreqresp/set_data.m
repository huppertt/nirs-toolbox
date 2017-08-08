function data = set_data(this, data)
%SET_DATA   PreSet function for the 'data' property.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

% Determine if the input is a matrix; if row make it column.
[data,nfft,nchans] = checkinputsigdim(data);

% When the data changes we have to assume that the stored spectrum
% information no longer applies.
this.Metadata.setsourcespectrum([]);

% [EOF]

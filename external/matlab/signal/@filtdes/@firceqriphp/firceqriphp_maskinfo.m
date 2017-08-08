function cmd = firceqriphp_maskinfo(hObj, d)
%FIRCEQRIPHP_MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = ft_maskinfo(hObj, d);

cmd.bands{1}.drawpatch    = false;
cmd.bands{1}.magfcn       = 'rolloff';
cmd.bands{1}.slope        = get(d, 'stopbandSlope');
cmd.bands{1}.drawfreqbars = false;

% [EOF]

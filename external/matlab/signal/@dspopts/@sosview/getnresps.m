function n = getnresps(this, Hd)
%NRESPS   Return the number of responses for the given filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

switch lower(this.View)
    case 'complete'
        n = 1;
    case {'cumulative', 'individual'}
        n = nsections(Hd);
    case 'userdefined'
        % Call trim custom so we throw away any responses that have all of
        % their indexes greater than nsections of Hd.
        n = length(trimcustom(this, Hd));
end

% [EOF]

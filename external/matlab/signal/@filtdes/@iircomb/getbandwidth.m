function varargout = getbandwidth(h)
%GETBANDWIDTH Returns the bandwidth

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

switch lower(h.TransitionMode)
    case 'bandwidth'
        if nargout == 1,
            varargout = {get(h, 'Bandwidth')};
        else
            varargout = {'BW', getmcode(h, 'Bandwidth'), '', ''};
        end
    case 'q'
        if nargout == 1,
            varargout = {2*getnyquist(h)/h.order/h.q};
        else
            [fs, fsstr] = getfsstr(h);
            varargout = {'Q', getmcode(h, 'Q'), 'Q-factor', ...
                    sprintf('\nBW = 2*%s/%d/Q;\n',fsstr, h.order)};
        end
end


% [EOF]

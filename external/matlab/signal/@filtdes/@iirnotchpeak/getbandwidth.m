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
        switch lower(h.ResponseType)
            case 'notching'
                if nargout == 1,
                    varargout = {h.fnotch/h.q};
                else
                    varargout = {'Q', getmcode(h, 'Q'), '', sprintf('\nBW = Fnotch/Q;\n')};
                end
            case 'peaking'
                if nargout == 1,
                    varargout = {h.fpeak/h.q};
                else
                    varargout = {'Q', getmcode(h, 'Q'), '', sprintf('\nBW = Fpeak/Q;\n')};
                end
        end
end

% [EOF]
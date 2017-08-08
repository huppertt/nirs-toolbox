function varargout = dispstr(this, varargin)
%DISPSTR   Convert the coefficients to the display

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Ignore the format, can only work fixed.
if ischar(varargin{end})
    fmt = varargin{end};
    varargin(end) = [];
else
    fmt = 'dec';
end

switch lower(fmt(1:3))
    case 'dec'
        for indx = 1:length(varargin)
            varargout{indx} = signal_num2str(varargin{indx});
        end
    case 'hex'
        for indx = 1:length(varargin)
            [rows, cols] = size(varargin{indx});
            if cols > 1 && rows > 1
                str = [];
                for jndx = 1:cols
                    str = [str repmat('  ', rows, 1) num2hex(varargin{indx}(:,jndx))];
                end
                str(:, 1:2) = [];
            else
                str = num2hex(varargin{indx});
            end
            varargout{indx} = str;
        end
    case 'bin'
        if isfixptinstalled
            q = quantizer(getarithmetic(this));
            for indx = 1:length(varargin)
                [rows, cols] = size(varargin{indx});
                if cols > 1 && rows > 1
                    str = [];
                    for jndx = 1:cols
                        str = [str repmat('  ', rows, 1) num2bin(q, varargin{indx}(:,jndx))];
                    end
                    str(:, 1:2) = [];
                else
                    str = num2bin(q, varargin{indx});
                end
                varargout{indx} = str;
            end
        else
            error(message('signal:dfilt:filterquantizer:dispstr:invalidFormat'));
        end
end

% [EOF]

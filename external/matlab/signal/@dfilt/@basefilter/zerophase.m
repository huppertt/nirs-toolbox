function varargout = zerophase(this, varargin)
%ZEROPHASE Zero-phase response of a discrete-time filter.
%   [Hr,W] = zerophase(Hd,N) returns length N vectors Hr and W containing
%   the zero-phase response of the discrete-time filter Hd, and the
%   frequencies (in radians) at which it is evaluated. The zero-phase
%   response is evaluated at N points equally spaced around the upper half
%   of the unit circle.  For an FIR filter where N is a power of two, the
%   computation is done faster using FFTs.  If you don't specify N, it
%   defaults to 8192.
%
%   [Hr,W] = ZEROPHASE(Hb) returns a matrix H if Hb is a vector.  Each
%   column of the matrix corresponds to each filter in the vector.  If a
%   row vector of frequency points is specified, each row of the matrix
%   corresponds to each filter in the vector.
%
%   For additional parameters, see SIGNAL/PHASEZ.
%
%   See also DFILT, SIGNAL/ZEROPHASE.

%   Author: V. Pellissier, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin > 1 && ischar(varargin{1}) && ~any(strcmpi(varargin{1}, {'whole', 'half'}))
    newobj = true;
else
    newobj = false;
end

if newobj
    hopts = uddpvparse('dspopts.freqresp', varargin{:});

    inputs = hopts.freqzinputs;

    [h,w,p,opts] = calculate_zerophase(this, inputs{:});

    opts = {};
    if ~hopts.NormalizedFrequency
        opts = {'Fs', hopts.Fs};
    end

    if strcmpi(hopts.FrequencySpecification, 'NFFT')
        opts = {opts{:}, 'SpectrumRange', hopts.SpectrumRange};
    end

    h = dspdata.zerophase(h, w, opts{:});

    if nargout
        varargout = {h, opts};
    else
        plot(h);
    end
elseif nargout,

    [h,w,p,opts] = calculate_zerophase(this, varargin{:});
    varargout = {h,w,p,opts};
else,
    [this, opts] = freqzparse(this, varargin{:});
    opts.magnitude = 'Zero-phase';
    opts.phase     = 'Continuous Phase';
    fvtool(this, opts);
end

% -------------------------------------------------------------------------
function [h, w, p, opts] = calculate_zerophase(this, varargin)

h = [];
p = [];

% Do not use BASE_RESP since zerophase needs 3 outputs, not 2.  It also
% requires the opts structure.  This is faster than using BASE_RESP
% with 2 different methods to overload as functional ZEROPHASE would
% have to be called twice.
for indx = 1:length(this),
    Hd = dispatch(this(indx));
    for jndx = 1:length(Hd)
        [ht, w, pt, opts] = computezerophase(Hd(jndx), varargin{:});

        if isempty(h),
            h = ht;
            p = pt;
        else

            % If the number of columns == 1, expand in columns
            if size(ht, 2) == 1,
                h(:, end+1) = ht;
                p(:, end+1) = pt;
            else
                h(end+1, :) = ht;
                p(end+1, :) = pt;
            end
        end
    end
end

% [EOF]

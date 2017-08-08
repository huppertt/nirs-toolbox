function varargout = freqresp(this, varargin)
%FREQRESP   Discrete-time filter frequency response.
%   H = FREQRESP(H) returns the complex frequency response object H.
%
%   H = FREQRESP(Hb, Param1, Value1, Param2, Value2, etc) returns the
%   response object given the options in P-V pairs.  Valid options are:
%   Parameter                Default     Description/Valid values
%   ---------                -------     ------------------------
%   'NormalizedFrequency'    true        
%   'Fs'                     1          Not used when NormalizedFrequency
%                                       is true.
%   'FrequencySpecification' 'NFFT'     {'NFFT', 'FrequencyVector'}
%   'NFFT'                   8192       Not used when FrequencySpecification
%                                       is set to 'FrequencyVector'
%   'FrequencyVector'       linspace(0, 1, 512)
%   'SpectrumRange'         'Half'      {'Whole', 'Half'}
%   'CenterDC'              false
%
%   These options are all contained in the DSPOPTS.FREQRESP object.
%
%   See also DFILT, SIGNAL/FREQZ.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

hopts = uddpvparse('dspopts.freqresp', varargin{:});

inputs = freqzinputs(hopts);

[h,w] = base_resp(this, 'computefreqz', inputs{:});

opts = {};
if ~hopts.NormalizedFrequency
    opts = {'Fs', hopts.Fs};
end

if strcmpi(hopts.FrequencySpecification, 'NFFT')
    opts = {opts{:}, 'SpectrumRange', hopts.SpectrumRange};
end

h = dspdata.freqresp(h, w, opts{:});

if nargout,
    varargout = {h};
else
    plot(h);
end

% [EOF]

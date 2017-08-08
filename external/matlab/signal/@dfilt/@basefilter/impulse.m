function varargout = impulse(this, varargin)
%IMPULSE Impulse response of digital filter
%   H = IMPULSE(Hb) returns the impulse response object H.
%
%   H = IMPULSE(Hb, PV Pairs) returns the response object H based on PV
%   Pairs.  Valid options are:
%   Parameter                Default     Description/Valid values
%   ---------                -------     ------------------------
%   'NormalizedFrequency'    true        
%   'Fs'                     1          Not used when NormalizedFrequency
%                                       is true.
%   'LengthOption'           'Default'  {'Default', 'Specified'}
%   'Length'                 20         Not used when LengthOption is set
%                                       to 'Default'.
%
%   These options are all contained in the DSPOPTS.TIMERESP object.
%
%   For additional parameters, see SIGNAL/IMPZ.
%
%   See also DFILT, SIGNAL/IMPZ.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

feature('TimeSeriesTools',1);

hopts = uddpvparse('dspopts.timeresp', varargin{:});

inputs = oldinputs(hopts);

[y, t] = base_resp(this, 'computeimpz', inputs{:});

h = tsdata.timeseries(y, t);

if nargout,
    varargout = {h};
else
    hax = newplot;
    
    t = get(h, 'Time');
    
    if hopts.NormalizedFrequency
        xunits = getString(message('signal:dfilt:dfilt:Samples'));
    else
        [t, m, xunits] = engunits(t, 'time', 'latex'); %#ok
        xunits = [xunits 's'];
    end
    stem(hax, t, h.Data);
    xlabel(hax, sprintf('%s',getString(message('signal:dfilt:dfilt:Time', xunits))));
    ylabel(hax, getString(message('signal:dfilt:dfilt:Amplitude')));
    title(hax, getString(message('signal:dfilt:dfilt:ImpulseResponse')));
%     hs = stem(h);
%     title(ancestor(hs, 'axes'), getString(message('signal:dfilt:dfilt:ImpulseResponse')));
end

% [EOF]

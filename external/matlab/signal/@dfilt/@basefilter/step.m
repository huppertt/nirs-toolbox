function varargout = step(this, varargin)
%STEP   Step response.
%   H = STEP(Hb) computes the step response object H.
%
%   For additional parameters see DFILT.BASEFILTER/IMPULSE.
%
%   See also DFILT, SIGNAL/STEPZ.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

feature('TimeSeriesTools',1);

hopts = uddpvparse('dspopts.timeresp', varargin{:});

inputs = oldinputs(hopts);

[y, t] = base_resp(this, 'computestepz', inputs{:});

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
    title(hax, getString(message('signal:dfilt:dfilt:StepResponse')));

%     hs = stem(h);
%     title(ancestor(hs, 'axes'), getString(message('signal:dfilt:dfilt:StepResponse')));
end

% [EOF]

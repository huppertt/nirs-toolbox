function varargout = plot(this)
%PLOT   Plot the response.

%   Copyright 1988-2014 The MathWorks, Inc.

if length(this) > 1
    error(message('signal:dspdata:abstractfreqresp:plot:InvalidInputs'));
end

normfreq = get(this, 'NormalizedFrequency');

% Determine the frequency range to plot.
freqrange = 'whole';
if ishalfnyqinterval(this)
    freqrange = 'half';
end
centerdc = getcenterdc(this);

% Create a new plot or reuse an available one.
hax = newplot;
    
% Get the data from this object.
[H, W] = getdata(this,isdensity(this),plotindb(this),normfreq,freqrange,centerdc);

% Set up the xlabel.
if normfreq
    W    = W/pi;
    xlbl = getfreqlbl('rad/sample');
else
    [W, ~, xunits] = engunits(W);
    xlbl = getfreqlbl([xunits 'Hz']);
end

% Plot the data.
h = line(W, H, 'Parent', hax);

if((strcmp(this.Name, 'Power Spectral Density') || strcmp(this.Name, 'Mean-Square Spectrum')) && ~isempty(this.ConfInterval))
    CI = this.ConfInterval;
    CL = this.ConfLevel;   
    Hc = db(CI,'power');
    
    % Plot the Confidence Intervals.
    for i=1:size(H,2)
        baseColor = get(h(i,1),'Color');
        h(i,2) = line(W, Hc(:,2*i-1),'Color',baseColor,'LineStyle','-.','Parent',hax);
        h(i,3) = line(W, Hc(:,2*i),'Color',baseColor,'LineStyle','-.','Parent', hax);
    end
    
    % convert to row vector for backwards compatibility
    h = h(:)';
    
    % re-order the children so first two legend entries are 'correct'.
    hc = get(hax,'Children');
    if numel(hc)==numel(h)
        set(hax,'Children',reshape(reshape(hc,numel(hc)/3,3)',1,numel(hc)));
    end

    % re-save as rows for backwards compatibility.
    h = h';
    
    if strcmp(this.Name, 'Power Spectral Density')
      Estimate = getString(message('signal:dspdata:abstractfreqresp:plot:PowerSpectralDensity'));
    else
      Estimate = getString(message('signal:dspdata:abstractfreqresp:plot:PowerSpectrum'));
    end
    Interval = getString(message('signal:dspdata:abstractfreqresp:plot:ConfidenceIntervalPct',num2str(CL*100)));
    legend(Estimate,Interval,'Location','best');    
end

xlabel(hax, xlbl);

% Set up the ylabel
ylabel(hax, getylabel(this));

title(hax, getTranslatedString('signal:dspdata:dspdata',gettitle(this)));

set(hax, 'Box', 'On', ...
    'XGrid', 'On', ...
    'YGrid', 'On', ...
    'XLim', [min(W) max(W)]);

% Ensure axes limits are properly cached for zoom/unzoom
resetplotview(hax,'SaveCurrentView'); 

if nargout
    varargout = {h};
end

% [EOF]

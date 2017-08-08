function timezplot(t,h,Fs,tstr)
%TIMEZPLOT Plot the time response data
%   TIMEZPLOT(T,H,FS,TITLE) Plot the time response represented by the time
%   vector T and the response vector H.  FS is the sampling frequency and
%   defaults to 1.  TITLE is the prefix to ' Response' for the axes title
%   and defaults to ''.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(2,4,nargin,'struct'));

if nargin < 3, Fs   = 1; end
if nargin < 4, tstr = ''; end

plotimag = any(~isreal(h));
if plotimag,
    subplot(2,1,1);
    ylbl = getString(message('signal:timezplot:Real'));
else
    ylbl = '';
end

if Fs == 1, xlbl = ['n (' getString(message('signal:timezplot:Samples')) ')'];
else,       xlbl = ['nT (' getString(message('signal:timezplot:Seconds')) ')']; end

lclstem(t, real(h), xlbl, ylbl);
title([tstr ' ' getString(message('signal:timezplot:Response'))]);

if plotimag,
    subplot(2,1,2);
    lclstem(t, imag(h), xlbl, getString(message('signal:timezplot:Imaginary')));
end

% Set the figure's next plot to replace to match the functionality of
% freqz.  This is necessary so subsequent plots to the figure replace
% both axes for the complex case.
set(gcf,'nextplot','replace');   


% -----------------------------------------------------------
function lclstem(t, h, xlbl, ylbl)

stem(t,h,'filled');

xlabel(xlbl); 
ylabel([ylbl ' ' getString(message('signal:timezplot:Amplitude'))]);

xlimits = [t(1) t(end)];
if ~isequal(xlimits(1), xlimits(2)),
    set(gca,'xlim',xlimits);
end


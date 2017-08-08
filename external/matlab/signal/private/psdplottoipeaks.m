function psdplottoipeaks(fscale, fundPow, fundFreq, imodPow, imodFreq)
%PSDPLOTTOIPEAKS Helper function for plotting power and psd estimates
%   used by TOI function

%   Copyright 2013 The MathWorks, Inc.

opts = {'EdgeColor','k','BackgroundColor','w', ...
        'VerticalAlignment','bottom','HorizontalAlignment','right'};

text(fundFreq(1)*fscale,fundPow(1),'F_{1}',opts{1:8},'FontSize',8);
text(fundFreq(2)*fscale,fundPow(2),'F_{2}',opts{1:6},'FontSize',8);
text(imodFreq(1)*fscale,imodPow(1),'2 F_{1} - F_{2}',opts{1:8},'FontSize',8);
text(imodFreq(2)*fscale,imodPow(2),'2 F_{2} - F_{1}',opts{1:6},'FontSize',8);
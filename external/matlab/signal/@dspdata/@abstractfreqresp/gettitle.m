function titlestr = gettitle(this)
%GETTITLE   Get the title.

%   Author(s): J. Schickler
%   Copyright 2004-2006 The MathWorks, Inc.

% Set up the title.
titlestr = getname(this); % Handle PSD, MSS, and Pseudospectrum

infoStruc = info(this);
if strcmpi(infoStruc.EstimationMethod, 'unknown')
    % Remove "Estimation" from the title, when spectrum objects didn't
    % create the DSPDATA object.
    idx = regexpi(titlestr,' Estimate','once');
    titlestr = titlestr(1:idx-1);
else
    titlestr = sprintf('%s %s',infoStruc.EstimationMethod, titlestr);
end

% [EOF]

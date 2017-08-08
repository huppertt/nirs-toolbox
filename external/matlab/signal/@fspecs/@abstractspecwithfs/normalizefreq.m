function normalizefreq(h,boolflag,Fs)
%NORMALIZEFREQ Normalize frequency specifications.

%   Copyright 2003-2011 The MathWorks, Inc.

% Check for error condition first
if (nargin > 2) && boolflag,
  error(message('signal:fspecs:abstractspecwithfs:normalizefreq:FsnotAllowed'));
end

cachecurrentnormalizedfreq(h)

% Check for early return condition next
if nargin < 2,
    boolflag = true;
end
if boolflag && h.NormalizedFrequency,
    return;
end

oldFs = h.privFs;
% Assign Fs if specified
if (nargin > 2) && ~boolflag,
    h.privFs = Fs;
end

oldnormfreq = h.NormalizedFrequency;
h.privNormalizedFreq = boolflag;

if (~h.privFs == 0),
    p = props2normalize(h);
    if boolflag,
        for n = 1:length(p),
            set(h,p{n},2*get(h,p{n})/h.privFs);
        end
    else
        if ~oldnormfreq
           % If normalized frequency was already false set Fs so that specs
           % are recomputed correctly
           cf = h.privFs/oldFs; % Correction factor
        else
            cf = h.privFs*0.5;
        end
        for n = 1:length(p),
            set(h,p{n},cf*get(h,p{n}));
        end
    end
end

normalizetime(h,oldFs,oldnormfreq);

% [EOF]

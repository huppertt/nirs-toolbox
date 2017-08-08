function [msg,nfft,Fs,window,noverlap,p,dflag,msgobj] = psdchk(P,x,y)
%PSDCHK Helper function for PSD, CSD, COHERE, and TFE.
%   [msg,nfft,Fs,window,noverlap,p,dflag]=PSDCHK(P,x,y) takes the cell 
%   array P and uses each element as an input argument.  Assumes P has 
%   between 0 and 7 elements which are the arguments to psd, csd, cohere
%   or tfe after the x (psd) or x and y (csd, cohere, tfe) arguments.
%   y is optional; if given, it is checked to match the size of x.
%   x must be a numeric vector.
%   Outputs:
%     msg - error message, [] if no error
%     nfft - fft length
%     Fs - sampling frequency
%     window - window vector
%     noverlap - overlap of sections, in samples
%     p - confidence interval, [] if none desired
%     dflag - detrending flag, 'linear' 'mean' or 'none'

%   Author(s): T. Krauss, 10-28-93
%   Copyright 1988-2002 The MathWorks, Inc.

msg = [];
msgobj = [];

if length(P) == 0 
% psd(x)
    nfft = min(length(x),256);
    window = hanning(nfft);
    noverlap = 0;
    Fs = 2;
    p = [];
    dflag = 'none';
elseif length(P) == 1
% psd(x,nfft)
% psd(x,dflag)
    if isempty(P{1}),   dflag = 'none'; nfft = min(length(x),256); 
    elseif isstr(P{1}), dflag = P{1};       nfft = min(length(x),256); 
    else              dflag = 'none'; nfft = P{1};   end
    Fs = 2;
    window = hanning(nfft);
    noverlap = 0;
    p = [];
elseif length(P) == 2
% psd(x,nfft,Fs)
% psd(x,nfft,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}),   dflag = 'none'; Fs = 2;
    elseif isstr(P{2}), dflag = P{2};       Fs = 2;
    else              dflag = 'none'; Fs = P{2}; end
    window = hanning(nfft);
    noverlap = 0;
    p = [];
elseif length(P) == 3
% psd(x,nfft,Fs,window)
% psd(x,nfft,Fs,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    if isstr(P{3})
        dflag = P{3};
        window = hanning(nfft);
    else
        dflag = 'none';
        window = P{3};
        if length(window) == 1, window = hanning(window); end
        if isempty(window), window = hanning(nfft); end
    end
    noverlap = 0;
    p = [];
elseif length(P) == 4
% psd(x,nfft,Fs,window,noverlap)
% psd(x,nfft,Fs,window,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    window = P{3};
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isstr(P{4})
        dflag = P{4};
        noverlap = 0;
    else
        dflag = 'none';
        if isempty(P{4}), noverlap = 0; else noverlap = P{4}; end
    end
    p = [];
elseif length(P) == 5
% psd(x,nfft,Fs,window,noverlap,p)
% psd(x,nfft,Fs,window,noverlap,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    window = P{3};
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isempty(P{4}), noverlap = 0; else noverlap = P{4}; end
    if isstr(P{5})
        dflag = P{5};
        p = [];
    else
        dflag = 'none';
        if isempty(P{5}), p = .95;    else    p = P{5}; end
    end
elseif length(P) == 6
% psd(x,nfft,Fs,window,noverlap,p,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    window = P{3};
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isempty(P{4}), noverlap = 0; else noverlap = P{4}; end
    if isempty(P{5}), p = .95;    else    p = P{5}; end
    if isstr(P{6})
        dflag = P{6};
    else
        msgobj = message('signal:psdchk:MustBeString','DFLAG');
        msg = getString(msgobj); 
        return
    end
end

% NOW do error checking
if (nfft<length(window)), 
    msgobj = message('signal:psdchk:WindowTooBig');
    msg = getString(msgobj);
end
if (noverlap >= length(window)),
    msgobj = message('signal:psdchk:OverlapTooBig');
    msg = getString(msgobj);
end
if (nfft ~= abs(round(nfft)))|(noverlap ~= abs(round(noverlap))),
    msgobj = message('signal:psdchk:NeedPosIntsForNFFTandNOVERLAP');
    msg = getString(msgobj);
end
if ~isempty(p),
    if (prod(size(p))>1)|(p(1,1)>1)|(p(1,1)<0),
        msgobj = message('signal:psdchk:NeedNormalizedConfidence');
        msg = getString(msgobj);
    end
end
if min(size(x))~=1 | ~isnumeric(x) | length(size(x))>2
    msgobj = message('signal:psdchk:RequiresVectorInput');
    msg = getString(msgobj);
end
if (nargin>2) & ( (min(size(y))~=1) | ~isnumeric(y) | length(size(y))>2 )
    msgobj = message('signal:psdchk:RequiresVectorInput');
    msg = getString(msgobj);
end
if (nargin>2) & (length(x)~=length(y))
    msgobj = message('signal:psdchk:MismatchedLength','X','Y');
    msg = getString(msgobj);
end

dflag = lower(dflag);
if strncmp(dflag,'none',1)
      dflag = 'none';
elseif strncmp(dflag,'linear',1)
      dflag = 'linear';
elseif strncmp(dflag,'mean',1)
      dflag = 'mean';
else
    msgobj = message('signal:psdchk:InvalidDFLAG');
    msg = getString(msgobj);
end

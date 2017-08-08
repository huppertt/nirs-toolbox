function this = differentiatormeas(hfilter, varargin)
%DIFFERENTIATORMEAS   Construct a DIFFERENTIATORMEAS object.

%   Author(s): J. Schickler
%   Copyright 2005-2006 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.differentiatormeas;

minfo = parseconstructorinputs(this, hfilter, varargin{:});

this.Fpass = get_diff_fpass(hfilter, minfo.Fpass, minfo.Apass);
this.Fstop = get_diff_fstop(hfilter, this.Fpass, minfo.Fstop, minfo.Astop);

% Apass represents the passband ripple.
this.Apass = get_diff_apass(hfilter, minfo.Fpass, minfo.Apass,this.Fs);

% Astop represents stopband attenuation.
this.Astop = get_diff_astop(hfilter, minfo.Fpass, minfo.Fstop, minfo.Astop,this.Fs);


%--------------------------------------------------------------------------
function measured_fpass = get_diff_fpass (hd, fpass, apass)

measured_fpass = [];

% return fpass if it is specified
if (~isempty(fpass))
    measured_fpass = fpass;
% These two cases are to support type IV Differentiators in which case we
% know the Fpass
elseif (~isempty(apass))
    measured_fpass = 1;
elseif ((isempty(apass) && isempty(fpass)))
    measured_fpass = 1;
end


%--------------------------------------------------------------------------
function measured_fstop = get_diff_fstop (hd, fpass, fstop, astop)

measured_fstop = [];

% Return fstop if it is specified; Will never be in a situation that Astop
% was specified when Fstop was not.
if (~isempty(fstop))
    measured_fstop = fstop;
end


%--------------------------------------------------------------------------
function measured_apass = get_diff_apass (hd, fpass, apass, Fs)

if isempty(fpass),
    fpass = 1;
elseif isnumeric(Fs),
    fpass = fpass/Fs*2;
end

measured_apass = [];
wpass = fpass*pi;

N = 4096;

% Find theoretical Apass
r = hd.getratechangefactors;
htheo = r(1)*linspace(0,wpass,N);

% For a differentiator, the max Apass will be a wpass
hact = freqz(hd, linspace(0, wpass, N));

measured_apass = max(db(htheo(2:end))-db(hact(2:end)))-min(db(htheo(2:end))-db(hact(2:end)));

%--------------------------------------------------------------------------
function measured_astop = get_diff_astop (hd, fpass, fstop, astop, Fs)

measured_astop = [];
if isnumeric(Fs),
    fpass = fpass/Fs*2;
    fstop = fstop/Fs*2;
end
wstop = fstop*pi;
wpass = fpass*pi;

% Not all fspecs have a FilterOrder property, we need to determine this so
% that we only compute Astop for type III (even order) filters
ord = order(hd);
isTypeIII = false;
if ~rem(ord,2)
    isTypeIII = true;
end

% Calculate astop if fstop and if type III only (since type IV doesn't have
% a stopband)
if (~isempty(fstop) && isTypeIII)
    % Calculate astop over stopband
    h_maxpass = find_max(hd, 0, wpass);
    h_maxstop = find_max(hd, wstop, pi);
    measured_astop = db(h_maxpass) - db(h_maxstop);

    % else return specified astop if available
elseif (~isempty(astop))
    measured_astop = astop;
end


%--------------------------------------------------------------------------
function y = find_max (hd, lo, hi)

N = 1024;
w_lo = lo;
w_hi = hi;

% calculate max of response from w_lo to w_hi
if (w_lo == 0) && (w_hi == pi)
    [h,w] = freqz(hd, N);
    
else
    [h,w] = freqz(hd, linspace(w_lo, w_hi, N));
end

[y,idx] = max(abs(h));
if ((idx == 1) || (idx == N))
    return;

else
    w_lo = w(max(1, idx-1));
    w_hi = w(min(N, idx+1));
end

% repeat for finer resolution
h = freqz(hd, linspace(w_lo, w_hi, N));
y = max(abs(h));

% [EOF]

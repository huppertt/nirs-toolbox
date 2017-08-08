function [h,a] = fircls1(n,wo,dp,ds,varargin)
%FIRCLS1 Low & high pass FIR filter design by constrained least-squares.
%   B = FIRCLS1(N,WO,DEVP,DEVS) is a length N+1 linear-phase lowpass FIR
%   filter with cut-off frequency W0 and maximum band deviations or ripples
%   (in linear units) DEVP and DEVS.
%
%   W0 is between 0 and 1 (1 corresponds to half the sampling frequency),
%   DEVP is the maximum passband deviation or ripple (in linear units) from
%   1, and DEVS is the maximum stopband deviation or ripple (in linear
%   units) from 0.
%
%   B = FIRCLS1(N,WO,DEVP,DEVS,'high') is a highpass filter.  For this
%   case, the order N must be even.  Otherwise, the order will be
%   incremented by one.
%
%   B = FIRCLS1(N,WO,DEVP,DEVS,WT) and B =
%   FIRCLS1(N,WO,DEVP,DEVS,WT,'high') meets a passband or stopband edge
%   requirement. If WT is in the passband, then the use of this parameter
%   ensures that |E(WT)| <= DEVP where E(w) is the error function.
%   Similarly, if WT is in the stopband, then the use of WT ensures that
%   |E(WT)| <= DEVS. Note that in the design of very narrow band filters
%   with small DEVP and DEVS, there may not exist a filter of the given
%   length that meets these specifications.
%
%   B = FIRCLS1(N,WO,DEVP,DEVS,WP,WS,K) weights the square error in the
%   passband K times greater than that in the stopband. WP is the passband
%   edge of the L2 weight function and WS is the stopband edge (WP < W0 <
%   WS).  Use trailing WT and 'high' arguments to meet a requirement or
%   design a highpass filter as in the case with no weighting function,
%   e.g.
%      FIRCLS1(N,WO,DEVP,DEVS,WP,WS,K,WT,'high').
%   In the 'high' pass filter case, you must have WS < W0 < WP.
%
%   MONITORING THE DESIGN
%      For a textual progress report on the iteration, use a trailing
%      'trace' argument, e.g. FIRCLS1(N,W0,DEVP,DEVS,...,'trace').  To
%      display plots of the design iteration, use a trailing 'plots'
%      argument, FIRCLS1(N,W0,DEVP,DEVS,...,'plots').  For both text and
%      plots, use 'both'.
%
%   EXAMPLES
%      n = 55;
%      wo = 0.3;
%      dp = 0.02; ds = 0.008;
%      h = fircls1(n,wo,dp,ds);           % no L2 weights
%      wp = 0.28; ws = 0.32;
%      K = 10;
%      h1 = fircls1(n,wo,dp,ds,wp,ws,K);  % L2 weight
%      fvtool(h,1,h1,1);
%
%   See also FIRCLS, FIRLS, FIRPM, DESIGNFILT.

%   Reference:
%   I. W. Selesnick, M. Lang and C. S. Burrus, Constrained Least Square
%   Design of FIR Filters Without Specified Transition Bands, IEEE Trans.
%   on Signal Processing, Aug 1996, Vol 44, Number 8, pp. 1879-1892.

%   Copyright 1988-2013 The MathWorks, Inc.     

%  subprograms : local_max.m, frefine.m

% --------- check input parameters ---------

if n < 3,
    error(message('signal:fircls1:invalidFilterOrder'));
end

LOW = 0;
HIGH = 1;
pass_type = LOW;
EDGE = 0;
WEIGHTS = 0;
TEXT_PF = 0;
PLOT_PF = 0;

if nargin >= 7 & ~ischar(varargin{3}) %#ok
    WEIGHTS = 1;
    % Cast to enforce Precision Rules
    wp = signal.internal.sigcasttofloat(varargin{1},'double','fircls1',...
      'WP','allownumeric');
    ws = signal.internal.sigcasttofloat(varargin{2},'double','fircls1',...
      'WS','allownumeric');
    K = signal.internal.sigcasttofloat(varargin{3},'double','fircls1','K',...
      'allownumeric');
    
    varargin = varargin(4:end);
end

if nargin < 4
    error(message('signal:fircls1:Nargchk'))
elseif prod([size(n), size(wo), size(dp), size(ds)]) > 1
    error(message('signal:fircls1:NeedScalar'))
elseif ( wo <= 0) || ( wo >= 1)
    error(message('signal:fircls1:InvalidW0Range'))
elseif (dp <= 0) || (ds <= 0)
    error(message('signal:fircls1:MustBePositive'))
elseif length(varargin)==1
    if ischar(varargin{1})
        switch lower(varargin{1}(1))
            case 'h'
                pass_type = HIGH;
            case 't'    % as in 'text'
                TEXT_PF = 1;
            case 'p'    % as in 'plots'
                PLOT_PF = 1;
            case 'b'    % as in 'both'
                TEXT_PF = 1; PLOT_PF = 1;
            otherwise
                error(message('signal:fircls1:InvalidParam'))
        end
    else
      % Cast to enforce Precision Rules
      wt = double(varargin{1});
        EDGE = 1;
    end
elseif length(varargin)==2
    if ischar(varargin{1})
        switch lower(varargin{1}(1))
            case 'h'
                pass_type = HIGH;
            otherwise
                error(message('signal:fircls1:InvalidParam'))
        end
        switch lower(varargin{2}(1))
            case 't'    % as in 'text'
                TEXT_PF = 1;
            case 'p'    % as in 'plots'
                PLOT_PF = 1;
            case 'b'    % as in 'both'
                TEXT_PF = 1; PLOT_PF = 1;
            otherwise
                error(message('signal:fircls1:InvalidParam'))
        end
    else
      % Cast to enforce Precision Rules
      wt = double(varargin{1});
        EDGE = 1;
        switch lower(varargin{2}(1))
            case 'h'
                pass_type = HIGH;
            case 't'    % as in 'text'
                TEXT_PF = 1;
            case 'p'    % as in 'plots'
                PLOT_PF = 1;
            case 'b'    % as in 'both'
                TEXT_PF = 1; PLOT_PF = 1;
            otherwise
                error(message('signal:fircls1:InvalidParam'))
        end
    end
elseif length(varargin)==3
  % Cast to enforce Precision Rules
  wt =  signal.internal.sigcasttofloat(varargin{1},'double','fircls1','',...
    'allownumeric');
    EDGE = 1;
    if ischar(varargin{2})
        switch lower(varargin{2}(1))
            case 'h'
                pass_type = HIGH;
            otherwise
                error(message('signal:fircls1:InvalidParam'))
        end
    end
    if ischar(varargin{3})
        switch lower(varargin{3}(1))
            case 't'    % as in 'text'
                TEXT_PF = 1;
            case 'p'    % as in 'plots'
                PLOT_PF = 1;
            case 'b'    % as in 'both'
                TEXT_PF = 1; PLOT_PF = 1;
            otherwise
                error(message('signal:fircls1:InvalidParam'))
        end
    end
end
% Cast to enforce Precision Rules
n = signal.internal.sigcasttofloat(n,'double','fircls1','N','allownumeric');
wo = signal.internal.sigcasttofloat(wo,'double','fircls1','WO',...
  'allownumeric');
dp = signal.internal.sigcasttofloat(dp,'double','fircls1','DEVP',...
  'allownumeric');
ds = signal.internal.sigcasttofloat(ds,'double','fircls1','DEVS',...
  'allownumeric');

PASS = 1;
STOP = 2;
if EDGE
    if isempty(wt),
        error(message('signal:fircls1:Empty'));
    elseif numel(wt) > 1
        error(message('signal:fircls1:MustBeScalar', 'WT'))
    elseif (wt <= 0) || (wt >= 1)
        error(message('signal:fircls1:InvalidRange', 'WT'))
    elseif wt == wo
        error(message('signal:fircls1:SignalErr'))
    elseif wt < wo
        if pass_type == LOW
            edge_type = PASS;
        else
            edge_type = STOP;
        end
    else
        if pass_type == LOW
            edge_type = STOP;
        else
            edge_type = PASS;
        end
    end
end

if WEIGHTS
    if pass_type == LOW
        if (wp > wo) || (ws < wo)
            error(message('signal:fircls1:InvalidLowpassRange'))
        end
    else
        if (ws > wo) || (wp < wo)
            error(message('signal:fircls1:InvalidHighpassRange'))
        end
    end
end

if pass_type == LOW
    up = [1+dp  ds]; d1 = dp;
    lo = [1-dp -ds]; d2 = ds;
    mag = [1 0];
else
    up = [ ds 1+dp]; d1 = ds;
    lo = [-ds 1-dp]; d2 = dp;
    mag = [0 1];
end


[n,msg1,msg2,msgobj] = firchk(n,1,pass_type);
if ~isempty(msg1), error(msgobj); end;
if ~isempty(msg2), warning(msgobj); end;

n = n+1;  % convert order to length for this routine
if rem(n,2) == 1
    parity = 1;
else
    parity = 0;
end

% --------- start algorithm ---------

wo = wo*pi;
if WEIGHTS
    wp = wp*pi;
    ws = ws*pi;
end
if EDGE, wt = wt*pi; end
L = 2^ceil(log2(5*n));
r = sqrt(2);
w = (0:L)*pi/L;                  % w includes both 0 and pi
q = round(wo*L/pi);
u = [up(1)*ones(q,1); up(2)*ones(L+1-q,1)];
l = [lo(1)*ones(q,1); lo(2)*ones(L+1-q,1)];
if parity  == 1
    m = (n-1)/2;
    if WEIGHTS
        c = 2*[wp/r; (sin(wp*(1:m))./(1:m))']/pi;
    else
        c = 2*[wo/r; (sin(wo*(1:m))./(1:m))']/pi;
    end
    if pass_type == HIGH
        c = -c;
        c(1) = c(1)+r;
    end
    Z = zeros(2*L-n,1);
    v = 0:m;
    tt = 1 - 1/r;
    NP = m+1;	% NP : number of parameters
else
    m = n/2;
    v = (1:m)-1/2;
    % c = [4*((sin(wo*[2*[1:m]-1]/2))./(2*[1:m]-1))/pi]';
    if WEIGHTS
        c = (2*sin(wp*v)./(v*pi))';
    else
        c = (2*sin(wo*v)./(v*pi))';
    end
    Z = zeros(4*L,1);
    tt = 0;
    NP = m;		% NP : number of parameters
end


if WEIGHTS
    % ----- construct R matrix --------------------
    if rem(n,2) == 1
        % odd length symmetric filter
        v1 = 1:m;
        v2 = m:2*m;
        if pass_type == LOW
            tp = [wp+K*(pi-ws), (sin(v1*wp)-K*sin(v1*ws))./v1]/pi;
            hk = ((sin(v2*wp)-K*sin(v2*ws))./v2)/pi;
        else % pass_type == HIGH
            tp = [(pi-wp)+K*ws, (-sin(v1*wp)+K*sin(v1*ws))./v1]/pi;
            hk = ((-sin(v2*wp)+K*sin(v2*ws))./v2)/pi;
        end
        R = toeplitz(tp,tp) + hankel(tp,hk);
        R(1,:) = R(1,:)/r;
        R(:,1) = R(:,1)/r;
        Ri = inv(R);
    else
        % even length symmetric filter
        v1 = 1:(m-1);
        tp = [(wp+K*(pi-ws)), (sin(v1*wp)-K*sin(v1*ws))./v1]/pi;
        v1 = 1:m;
        v2 = m:2*m-1;
        tp2 = ((sin(v1*wp)-K*sin(v1*ws))./v1)/pi;
        hk2 = ((sin(v2*wp)-K*sin(v2*ws))./v2)/pi;
        R = toeplitz(tp,tp) + hankel(tp2,hk2);
        Ri = inv(R);
    end

    a = Ri*c;       % best L2 cosine coefficients
else
    a = c;          % best L2 cosine coefficients
end
SN = 1e-8;              % Small Number
% -------- BEGIN LOOPING --------------
while 1

    % calculate H
    if parity == 1
        H = fft([a(1)*r; a(2:m+1); Z; a(m+1:-1:2)]);
        H = real(H(1:L+1)/2);
    else
        Z(2:2:2*m) = a;
        Z(4*L-2*m+2:2:4*L) = a(m:-1:1);
        H = fft(Z);
        H = real(H(1:L+1)/2);
    end

    % find local maxima and minima
    kmax = local_max(H);
    kmin = local_max(-H);
    % if filter length is even, then remove w=pi from constraint set
    if parity == 0
        n1 = length(kmax);
        if kmax(n1) == L+1, kmax(n1) = []; end
        n2 = length(kmin);
        if kmin(n2) == L+1, kmin(n2) = []; end
    end

    % refine frequencies
    cmax = (kmax-1)*pi/L;      cmin = (kmin-1)*pi/L;
    cmax = frefine(a,v,cmax);  cmin = frefine(a,v,cmin);

    % insert wt into constraint set if necessary
    if EDGE
        if pass_type == LOW
            if edge_type == PASS
                w_key = max(cmax(cmax<wo));
                if wt > w_key
                    kmin = [kmin; 1];
                    cmin = [cmin; wt];
                end
            else % edge_type == STOP
                w_key = min(cmin(cmin>wo));
                if wt < w_key
                    kmax = [kmax; L];
                    cmax = [cmax; wt];
                end
            end
        else % pass_type == HIGH
            if edge_type == PASS
                w_key = min(cmax(cmax>wo));
                if wt < w_key
                    kmin = [kmin; L];
                    cmin = [cmin; wt];
                end
            else % edge_type == STOP
                w_key = max(cmin(cmin<wo));
                if wt > w_key
                    kmax = [kmax; 1];
                    cmax = [cmax; wt];
                end
            end
        end
    end

    % evaluate H at refined frequencies
    Hmax = cos(cmax*v)*a - tt*a(1);
    Hmin = cos(cmin*v)*a - tt*a(1);

    % determine new constraint set
    v1   = Hmax > u(kmax)-100*SN;
    v2   = Hmin < l(kmin)+100*SN;
    kmax = kmax(v1); kmin = kmin(v2);
    cmax = cmax(v1); cmin = cmin(v2);
    Hmax = Hmax(v1); Hmin = Hmin(v2);
    n1   = length(kmax);
    n2   = length(kmin);

    % plot figures
    if PLOT_PF
        wv = [cmax; cmin];
        Hv = [Hmax; Hmin];
        subplot(311)
        plot(w/pi,H),
        axis([0 1 -.2 1.2])
        hold on,
        plot(wv/pi,Hv,'o'),
        hold off
        subplot(312)
        plot(w/pi,H-mag(1)),
        hold on,
        plot(wv/pi,Hv-mag(1),'o'),
        hold off
        axis([0 wo/pi -2*d1 2*d1])
        subplot(313)
        plot(w/pi,H-mag(2)),
        hold on,
        plot(wv/pi,Hv-mag(2),'o'),
        hold off
        axis([wo/pi 1 -2*d2 2*d2])
        pause(.5)
    end

    % remove a constraint set frequency if necessary (if otherwise overdetermined)
    if (n1+n2) > NP
        if parity == 1
            H0 =  a(1)/sqrt(2) + sum(a(2:m+1));
            Hpi = a(1)/sqrt(2) + sum(a(3:2:m+1)) - sum(a(2:2:m+1));
            dH0dw = -sum(a(2:m+1)'.*((1:m).^2));
            dHpidw = sum(a(2:2:m+1)'.*((1:2:m).^2)) - sum(a(3:2:m+1)'.*((2:2:m).^2));
            if dH0dw > 0, E0 = lo(1) - H0;
            else          E0 = H0 - up(1); end
            if dHpidw > 0, Epi = lo(2) - Hpi;
            else           Epi = Hpi - up(2); end
        else % parity == 0
            % when length is even, we know that
            %  we must remove w = 0;
            E0 = 0; Epi = 1;
        end
        if E0 > Epi
            % remove w = pi from constraint set
            [temp1, k1] = max(kmin);
            [temp2, k2] = max(kmax);
            if temp1 < temp2
                % pi is in kmax
                kmax(k2) = [];
                cmax(k2) = [];
                Hmax(k2) = [];
                n1 = n1 - 1;
            else
                % pi is in kmin
                kmin(k1) = [];
                cmin(k1) = [];
                Hmin(k1) = [];
                n2 = n2 - 1;
            end
        else
            % remove w = 0 from constraint set
            [temp1, k1] = min(kmin);
            [temp2, k2] = min(kmax);
            if temp1 < temp2
                % 0 is in kmin
                kmin(k1) = [];
                cmin(k1) = [];
                Hmin(k1) = [];
                n2 = n2 - 1;
            else
                % 0 is in kmax
                kmax(k2) = [];
                cmax(k2) = [];
                Hmax(k2) = [];
                n1 = n1 - 1;
            end
        end
    end

    % check stopping criterion
    E  = max([Hmax-u(kmax); l(kmin)-Hmin; 0]);
    if TEXT_PF
        fprintf(1,'    Bound Violation = %15.13f  \n',E);
    end
    if E < SN
        break
    end

    % calculate new Lagrange multipliers
    if parity == 1
        G = [ones(n1,m+1); -ones(n2,m+1)].*cos([cmax; cmin]*(0:m));
        G(:,1) = G(:,1)/r;
    else
        G = [ones(n1,m); -ones(n2,m)].*cos([cmax; cmin]*((1:m)-1/2));
    end
    d = [u(kmax); -l(kmin)];

    if WEIGHTS
        mu = (G*Ri*G')\(G*Ri*c-d);
    else
        mu = (G*G')\(G*c-d);
    end

    % iteratively remove negative multiplier
    [min_mu,K] = min(mu);
    while min_mu < 0
        G(K,:) = [];
        d(K) = [];
        if WEIGHTS
            mu = (G*Ri*G')\(G*Ri*c-d);
        else
            mu = (G*G')\(G*c-d);
        end
        [min_mu,K] = min(mu);
    end

    % determine new cosine coefficients
    if WEIGHTS
        a = Ri*(c-G'*mu);
    else
        a = c-G'*mu;
    end

end

if parity == 1
    h = [a(m+1:-1:2)/2; a(1)/r; a(2:m+1)/2]';
else
    h = [a(m:-1:1); a]'/2;
end

if nargout > 1
    a = 1;
end

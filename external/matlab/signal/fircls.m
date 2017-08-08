function [h,a] = fircls(n,be,mag,up,lo,PF)
%FIRCLS  Linear-phase FIR filter design by constrained least-squares.
%   B = FIRCLS(N,F,A,DEV_UP,DEV_LO) is a length N+1 linear phase FIR filter
%   which approximates the desired piecewise constant frequency response in
%   F and A subject to the constraints given in DEV_UP and DEV_LO.
%
%   A is a vector describing the piecewise constant desired amplitude of
%   the frequency response. The length of A is the number of different
%   bands.
%
%   DEV_UP and DEV_LO are vectors of the same length of A. They give the
%   upper and lower maximum deviation or ripple (in linear units) for the
%   zerophase response in each band.
%
%   F is a vector of transition frequencies. These frequencies must be in
%   increasing order, must start with 0.0, and must end with 1.0. The
%   length of F should exceed the length of A by exactly 1.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass and
%   bandstop filters, N must be even.  Otherwise, N will be incremented by
%   one.
%
%   NOTE
%      By setting DEV_LO equal to 0 in the stopbands, a nonnegative
%      frequency response amplitude can be obtained. Such filters can be
%      spectrally factored to obtain minimum phase filters.
%
%   MONITORING THE DESIGN
%      For a textual progress report on the iteration, use a trailing
%      'trace' argument, e.g. FIRCLS(N,F,A,DEV_UP,DEV_LO,'trace').  To
%      display plots of the design iteration, use a trailing 'plots'
%      argument, FIRCLS(N,F,A,DEV_UP,DEV_LO,'plots').  For both text and
%      plots, use 'both'.
%
%   EXAMPLE
%      n   = 51;
%      f   = [0 0.4 0.8 1];
%      a = [0 1 0];
%      up  = [ 0.02 1.02  0.01];
%      lo  = [-0.02 0.98 -0.01];
%      b   = fircls(n,f,a,up,lo);
%
%   See also FIRCLS1, FIRLS, FIRPM, DESIGNFILT.

%   Reference:
%   I. W. Selesnick, M. Lang and C. S. Burrus, Constrained Least Square
%   Design of FIR Filters Without Specified Transition Bands, IEEE Trans.
%   on Signal Processing, Aug 1996, Vol 44, Number 8, pp. 1879-1892.

%   Copyright 1988-2013 The MathWorks, Inc.     

%   subprograms : local_max.m, frefine.m

% -------- check input parameters -------------

narginchk(5,6);

n = signal.internal.sigcasttofloat(n,'double','fircls','N',...
  'allownumeric');

if n < 3,
    error(message('signal:fircls:InvalidFilterOrder'));
end

% Cast to enforce precision rules
be = signal.internal.sigcasttofloat(be,'double','fircls','F',...
  'allownumeric');
mag = signal.internal.sigcasttofloat(mag,'double','fircls','A',...
  'allownumeric');
up = signal.internal.sigcasttofloat(up,'double','fircls','DEV_UP',...
  'allownumeric');
lo = signal.internal.sigcasttofloat(lo,'double','fircls','DEV_LO',...
  'allownumeric');

[n,msg1,msg2,msgobj] = firchk(n,be(end),mag);
if ~isempty(msg1), error(msgobj); end;
if ~isempty(msg2), warning(msgobj); end;

lm = length(mag);
lu = length(up);
ll = length(lo);
if (lm ~= lu) || (lm ~= ll) || (lu ~= ll)
    error(message('signal:fircls:UnequalLengths'))
elseif length(be) ~= (lm+1)
    error(message('signal:fircls:InvalidFDimensions'))
elseif be(1) ~= 0
    error(message('signal:fircls:InvalidFreqVecStart'))
elseif be(lm+1) ~= 1
    error(message('signal:fircls:InvalidFreqVecEnd'))
elseif any(diff(be)<=0)
    error(message('signal:fircls:InvalidFreqVecMonotonic'))
elseif any(up <= lo)
    error(message('signal:fircls:InvalidRangeUpGtLo'))
elseif any(lo > mag)
    error(message('signal:fircls:InvalidRangeLoLtA'))
elseif any(up < mag)
    error(message('signal:fircls:InvalidRangeUpGtA'))
end


n = n+1;  % convert order to length

TEXT_PF = 0;
PLOT_PF = 0;
if nargin == 6
    PF = lower(PF);
    switch PF(1)
        case 't'    % as in 'text'
            TEXT_PF = 1;
        case 'p'    % as in 'plots'
            PLOT_PF = 1;
        case 'b'    % as in 'both'
            TEXT_PF = 1; PLOT_PF = 1;
        otherwise,
            error(message('signal:fircls:InvalidParam'))
    end
end

L = 2^ceil(log2(3*n));
num_bands = length(mag);

if rem(n,2) == 0
    parity = 0;
    m = n/2;
else
    parity = 1;
    m = (n-1)/2;
end

r = sqrt(2);
be = be * pi;
% ----- calculate Fourier coefficients and upper ---
% ----- and lower bound functions ------------------
if parity == 1
    c = zeros(m+1,1);
    Z = zeros(2*L-1-2*m,1);
    v = 0:m;
    tt = 1-1/r;
else
    c = zeros(m,1);
    Z = zeros(4*L,1);
    v = (1:m)-1/2;
    tt = 0;
end
u = [];
l = [];
for k = 1:num_bands
    if parity == 1
        c = c + mag(k) * ...
            (2*[be(k+1)/r; [sin(be(k+1)*(1:m))./(1:m)]']/pi - ...
            2*[be(k)/r;   [sin(be(k)*(1:m))./(1:m)]']/pi); %#ok
    else
        c = c + mag(k) * ...
            ( [4*((sin(be(k+1)*[2*(1:m)-1]/2))./(2*(1:m)-1))/pi]' - ...
            [4*((sin(be(k) * [2*(1:m)-1]/2))./(2*(1:m)-1))/pi]' ) ;%#ok

    end
    q = round(L*(be(k+1)-be(k))/pi);
    u = [u; up(k)*ones(q,1)];
    l = [l; lo(k)*ones(q,1)];
end
flen = length(u);
if flen < L+1
    ov = ones(L+1-flen,1);
    u(flen+1:L+1) = up(num_bands)*ov;
    l(flen+1:L+1) = lo(num_bands)*ov;
elseif flen > L+1
    u = u(1:L+1);
    l = l(1:L+1);
end

w = (0:L)'*pi/L;
a = c;       % best L2 cosine coefficients
SN = 1e-8;   % Small Number
it = 0;
kmax = zeros(0,1); kmin = zeros(0,1);
cmax = zeros(0,1); cmin = zeros(0,1);
okmax = []; okmin = []; uvo = 0;
ocmax = []; ocmin = []; lvo = 0;

while 1

    if (uvo < -SN/10) || (lvo < -SN/10)
        % ----- include old extremal ----------------
        if uvo < lvo
            kmax = [kmax; okmax(k1)]; okmax(k1) = [];
            cmax = [cmax; ocmax(k1)]; ocmax(k1) = [];
        else
            kmin = [kmin; okmin(k2)]; okmin(k2) = [];
            cmin = [cmin; ocmin(k2)]; ocmin(k2) = [];
        end
    else
        % ----- calculate A -------------------------
        if parity == 1
            A = fft([a(1)*r;a(2:m+1);Z;a(m+1:-1:2)]);
            A = real(A(1:L+1))/2;
        else
            Z(2:2:2*m) = a;
            Z(4*L-2*m+2:2:4*L) = a(m:-1:1);
            A = fft(Z);
            A = real(A(1:L+1)/2);
        end

        if PLOT_PF
            cc = [cmax; cmin];
            occ = [ocmax; ocmin];
            if ~isempty(cc)
                Acc = cos(cc*v)*a - tt*a(1);
            else
                Acc = [];
            end
            if ~isempty(occ)
                Aocc = cos(occ*v)*a - tt*a(1);
            else
                Aocc = [];
            end

            subplot(length(be),1,1)
            plot(w/pi,A), hold on
            plot(cc/pi,Acc,'o')
            plot(occ/pi,Aocc,'x')
            hold off

            for k = 1:length(be)-1
                subplot(length(be),1,k+1);
                plot(w/pi,A), hold on
                plot(cc/pi,Acc,'o')
                plot(occ/pi,Aocc,'x')
                hold off
                axis([be(k)/pi be(k+1)/pi ...
                    mag(k)-1.5*(mag(k)-lo(k)) mag(k)+1.5*(up(k)-mag(k))])
                ylabel(sprintf('Band #%g',k))
            end
            xlabel('Frequency')
            drawnow
        end

        % ----- find extremals ----------------------
        okmax = kmax;              okmin = kmin;
        ocmax = cmax;              ocmin = cmin;
        kmax = local_max(A);       kmin = local_max(-A);
        if parity == 0
            n1 = length(kmax);
            if kmax(n1) == L+1, kmax(n1) = []; end
            n2 = length(kmin);
            if kmin(n2) == L+1, kmin(n2) = []; end
        end
        cmax = (kmax-1)*pi/L;      cmin = (kmin-1)*pi/L;
        cmax = frefine(a,v,cmax);  cmin = frefine(a,v,cmin);

        Amax = cos(cmax*v)*a - tt*a(1);
        Amin = cos(cmin*v)*a - tt*a(1);
        v1 = Amax > u(kmax)-SN;
        v2 = Amin < l(kmin)+SN;
        kmax = kmax(v1); kmin = kmin(v2);
        cmax = cmax(v1); cmin = cmin(v2);
        Amax = Amax(v1); Amin = Amin(v2);

        % ----- check stopping criterion ------------
        Eup = Amax-u(kmax); Elo = l(kmin)-Amin;
        E = max([Eup; Elo; 0]);
        if TEXT_PF
            fprintf(1,'  Bound Violation = %15.13f  \n',E);
        end
        if E < SN, break, end
    end

    % ----- calculate new multipliers -----------
    n1 = length(kmax);         n2 = length(kmin);
    if parity == 1
        G = [ones(n1,m+1); -ones(n2,m+1)] .* cos([cmax;cmin]*(0:m));
        G(:,1) = G(:,1)/r;
    else
        G = [ones(n1,m); -ones(n2,m)] .* cos([cmax;cmin]*v);
    end
    d  = [u(kmax); -l(kmin)];
    mu = (G*G')\(G*c-d);

    % ----- remove negative multiplier ----------
    [min_mu,K] = min(mu);
    while min_mu < 0
        G(K,:) = []; d(K) = [];
        mu = (G*G')\(G*c-d);
        if K > n1
            kmin(K-n1) = [];
            cmin(K-n1) = [];
            %if length(Amin)>=K-n1
            %   Amin(K-n1) = [];
            %end
            n2 = n2 - 1;
        else
            kmax(K) = [];
            cmax(K) = [];
            %Amax(K) = [];
            n1 = n1 - 1;
        end
        [min_mu,K] = min(mu);
    end

    % ----- determine new coefficients ----------
    a = c-G'*mu;

    if ~isempty(ocmax)
        oAmax = cos(ocmax*v)*a - tt*a(1);
        [uvo,k1] = min(u(okmax) - oAmax);
    else
        uvo = 0;
    end
    if ~isempty(ocmin)
        oAmin = cos(ocmin*v)*a - tt*a(1);
        [lvo,k2] = min(oAmin - l(okmin));
    else
        lvo = 0;
    end

    it = it + 1;
end

if parity == 1
    h = [a(m+1:-1:2); a(1)*r; a(2:m+1)]/2;
else
    h = [a(m:-1:1); a]/2;
end

h = h(:)';

if nargout > 1
    a = 1;
end

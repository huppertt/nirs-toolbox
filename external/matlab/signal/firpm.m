function [h,err,res] = firpm(order, ff, aa, varargin)
%FIRPM Parks-McClellan optimal equiripple FIR filter design.
%   B=FIRPM(N,F,A) returns a length N+1 linear phase (real, symmetric
%   coefficients) FIR filter which has the best approximation to the
%   desired frequency response described by F and A in the minimax sense. 
%   F is a vector of frequency band edges in pairs, in ascending order
%   between 0 and 1. 1 corresponds to the Nyquist frequency or half the
%   sampling frequency. At least one frequency band must have a non-zero
%   width. A is a real vector the same size as F which specifies the
%   desired amplitude of the frequency response of the resultant filter B.
%
%   The desired response is the line connecting the points (F(k),A(k)) and
%   (F(k+1),A(k+1)) for odd k; FIRPM treats the bands between F(k+1) and
%   F(k+2) for odd k as "transition bands" or "don't care" regions. Thus
%   the desired amplitude is piecewise linear with transition bands. The
%   maximum error is minimized.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass
%   and bandstop filters, N must be even.  Otherwise, N will be
%   incremented by one. Alternatively, you can use a trailing 'h' flag to
%   design a type 4 linear phase filter and avoid incrementing N.
%
%   B=FIRPM(N,F,A,W) uses the weights in W to weight the error. W has one
%   entry per band (so it is half the length of F and A) which tells
%   FIRPM how much emphasis to put on minimizing the error in each band
%   relative to the other bands.
%
%   B=FIRPM(N,F,A,'Hilbert') and B=FIRPM(N,F,A,W,'Hilbert') design filters
%   that have odd symmetry, that is, B(k) = -B(N+2-k) for k = 1, ..., N+1.
%   A special case is a Hilbert transformer which has an approx. amplitude
%   of 1 across the entire band, e.g. B=FIRPM(30,[.1 .9],[1 1],'Hilbert').
%
%   B=FIRPM(N,F,A,'differentiator') and B=FIRPM(N,F,A,W,'differentiator')
%   also design filters with odd symmetry, but with a special weighting
%   scheme for non-zero amplitude bands. The weight is assumed to be equal
%   to the inverse of frequency times the weight W. Thus the filter has a
%   much better fit at low frequency than at high frequency. This designs
%   FIR differentiators.
%
%   B=FIRPM(...,{LGRID}), where {LGRID} is a one-by-one cell array
%   containing an integer, controls the density of the frequency grid. The
%   frequency grid size is roughly LGRID*N/2*BW, where BW is the fraction
%   of the total band interval [0,1] covered by F. LGRID should be no less
%   than its default of 16. Increasing LGRID often results in filters which
%   are more exactly equiripple, at the expense of taking longer to
%   compute.
%
%   [B,ERR]=FIRPM(...) returns the maximum ripple height ERR.
%
%   [B,ERR,RES]=FIRPM(...) returns a structure RES of optional results
%   computed by FIRPM, and contains the following fields:
%
%      RES.fgrid: vector containing the frequency grid used in
%                 the filter design optimization
%        RES.des: desired response on fgrid
%         RES.wt: weights on fgrid
%          RES.H: actual frequency response on the grid
%      RES.error: error at each point on the frequency grid (desired - actual)
%      RES.iextr: vector of indices into fgrid of extremal frequencies
%      RES.fextr: vector of extremal frequencies
%
%   FIRPM is now a "function function", similar to CFIRPM, allowing you
%   to write a function which defines the desired frequency response.
%
%   B=FIRPM(N,F,@fresp,W) returns a length N+1 FIR filter which has the
%   best approximation to the desired frequency response as returned by the
%   function handle @fresp.  The function is called from within FIRPM using
%   the syntax:
%                    [DH,DW] = fresp(N,F,GF,W);
%   where:
%   N is the filter order.
%   F is the vector of frequency band edges which must appear monotonically
%     between 0 and +1, where 1 is the Nyquist  frequency.  The frequency
%     bands span F(k) to F(k+1) for k odd; the intervals  F(k+1) to F(k+2)
%     for k odd are "transition bands" or "don't care" regions during 
%     optimization.
%   GF is a vector of grid points which have been linearly interpolated
%     over each specified frequency band by FIRPM, and determines the
%     frequency grid at which the response function will be evaluated.
%   W is a vector of real, positive weights, one per band, for use
%     during optimization.  W is optional; if not specified, it is set
%     to unity weighting before being passed to 'fresp'.
%   DH and DW are the desired complex frequency response and
%     optimization weight vectors, respectively, evaluated at each
%     frequency in grid GF.
%
%   The predefined frequency response function handle for FIRPM is
%   @firpmfrf, but you can write your own.  See the help for
%   PRIVATE/FIRPMFRF for more information.
%
%   B=FIRPM(N,F,{@fresp,P1,P2,...},W) specifies optional arguments
%   P1, P2, etc., to be passed to the response function handle @fresp.
%
%   B=FIRPM(N,F,A,W) is a synonym for B=FIRPM(N,F,{@firpmfrf,A},W),
%   where A is a vector of response amplitudes at each band edge in F.
%
%   FIRPM normally designs symmetric (even) FIR filters. B=FIRPM(...,'h')
%   and B=FIRPM(...,'d') design antisymmetric (odd) filters. Each frequency
%   response function handle @fresp can tell FIRPM to design either an even
%   or odd filter in the absence of the 'h' or 'd' flags.  This is done
%   with:
%         SYM = fresp('defaults',{N,F,[],W,P1,P2,...})
%   FIRPM expects @fresp to return SYM = 'even' or SYM = 'odd'. If @fresp
%   does not support this call, FIRPM assumes 'even' symmetry.
%
%   % Example of a length 31 lowpass filter:
%   	h=firpm(30,[0 .1 .2 .5]*2,[1 1 0 0]);
%
%   % Example of a low-pass differentiator:
%   	h=firpm(44,[0 .3 .4 1],[0 .2 0 0],'differentiator');
%
%   % Example of a type 4 highpass filter:
%       h=firpm(25,[0 .4 .5 1],[0 0 1 1],'h');
%
%   See also FIRPMORD, CFIRPM, FIRLS, FIR1, FIR2, BUTTER, CHEBY1, CHEBY2,
%            ELLIP, FREQZ, FILTER, DESIGNFILT.

%   Author(s): L. Shure, 3-27-87
%          L. Shure, 6-8-88, revised
%          T. Krauss, 3-17-93, fixed hilbert bug 
%          T. Krauss, 3-2-97, consolidated grid generation, function-function

%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] "Programs for Digital Signal Processing", IEEE Press
%          John Wiley & Sons, 1979, pg. 5.1-1.
%     [2] "Selected Papers in Digital Signal Processing, II",
%          IEEE Press, 1976, pg. 97.

% Cast to enforce Precision rules
order = signal.internal.sigcasttofloat(order,'double','firpm','N','allownumeric');
ff = signal.internal.sigcasttofloat(ff,'double','firpm','F','allownumeric');

[nfilt,ff,grid,des,wt,ftype,sign_val,hilbert,neg] = firpminit(order, ff, aa, varargin{:});
% cast to enforce precision rules
wt = double(wt);
des = double(des);
grid = double(grid);

% Call actual design algorithm
[h,err,iext] = remezm(nfilt,ff/2,grid/2,des,wt,neg);

err = abs(err);
h = h(:).';  % make it a row
h = [h sign(.5-(neg))*h(length(h)-rem(nfilt,2):-1:1)];
h = h(length(h):-1:1);
if neg && ~hilbert
    h = -h;  %make sure differentiator has correct sign
end


%
% arrange 'results' structure
%
if nargout > 2 
    res.fgrid = grid(:);
    res.des = des(:);
    res.wt = wt(:);
    res.H = freqz(h,1,res.fgrid*pi);
    if neg  % asymmetric impulse response
        linphase = exp(sqrt(-1)*(res.fgrid*pi*(order/2) - pi/2));
    else
        linphase = exp(sqrt(-1)*res.fgrid*pi*(order/2));
    end
    if hilbert
        res.error = real(des(:) + res.H.*linphase);
    else
        res.error = real(des(:) - res.H.*linphase);
    end
    res.iextr = iext(1:end-1);
    res.fextr = grid(res.iextr);  % extremal frequencies
    res.fextr = res.fextr(:);
end

%%
function [h,dev,iext]  = remezm(nfilt,edge,grid,des,wt,neg)
% remezm function
% Inputs
%     nfilt - filter length
%     edge - vector of band edges (between 0 and .5)
%     grid - frequency grid (between 0 and .5)
%     des - desired function on frequency grid
%     wt - weight function on frequency grid
%     neg == 1 ==> antisymmetric imp resp,
%         == 0 ==> symmetric imp resp
% Outputs
%     h - coefficients of basis functions
%     dev - computed error
%     iext - indices of extremal frequencies

nbands = length(edge)/2;
jb = 2*nbands;
nodd = rem(nfilt,2);      % nodd == 1 ==> filter length is odd
% nodd == 0 ==> filter length is even
nfcns = fix(nfilt/2);
if nodd == 1 && neg == 0
    nfcns = nfcns + 1;
end

ngrid = length(grid);

if neg <= 0
    if nodd ~= 1
        des = des./cos(pi*grid);
        wt = wt.*cos(pi*grid);
    end
elseif nodd ~= 1
    des = des./sin(pi*grid);
    wt = wt.*sin(pi*grid);
else
    des = des./sin(2*pi*grid);
    wt = wt.*sin(2*pi*grid);
end
temp = (ngrid-1)/nfcns;
j=1:nfcns;
iext = fix([temp*(j-1)+1 ngrid])';
nz = nfcns + 1;

% Remez exchange loop
comp = [];
itrmax = 250;
devl = -1;
nzz = nz + 1;
niter = 0;
jchnge = 1;
jet = fix((nfcns - 1)/15) + 1;
ad = zeros(1,nz);
while jchnge > 0
    iext(nzz) = ngrid + 1;
    niter = niter + 1;
    if niter > itrmax
        break;
    end
    l = iext(1:nz)';
    x = cos(2*pi*grid(l));
    for nn = 1:nz
        ad(nn) = remezdd(nn, nz, jet, x);
    end
    add = ones(size(ad));
    add(2:2:nz) = -add(2:2:nz);
    dnum = ad*des(l)';
    dden = add*(ad./wt(l))';
    dev = dnum/dden;
    nu = 1;
    if dev > 0
        nu = -1;
    end
    dev = -nu*dev;
    y = des(l) + nu*dev*add./wt(l);
    if dev <= devl
        warning(message('signal:firpm:DidNotConverge',niter))
        break;
    end
    devl = dev;
    jchnge = 0;
    k1 = iext(1);
    knz = iext(nz);
    klow = 0;
    nut = -nu;
    j = 1;
    flag34 = 1;
    while j < nzz
        kup = iext(j+1);
        l = iext(j) + 1;
        nut = -nut;
        if j == 2
            y1 = comp;
        end
        comp = dev;
        flag = 1;
        if l < kup
            % gee
            c = ad./(cos(2*pi*grid(l))-x);
            err = (c*y'/sum(c) - des(l))*wt(l);
            dtemp = nut*err - comp;
            if dtemp > 0
                comp = nut*err;
                l = l + 1;
                while l < kup
                    % gee
                    c = ad./(cos(2*pi*grid(l))-x);
                    err = (c*y'/sum(c) - des(l))*wt(l);
                    dtemp = nut*err - comp;
                    if dtemp > 0
                        comp = nut*err;
                        l = l + 1;
                    else
                        break;
                    end
                end
                iext(j) = l - 1;
                j = j + 1;
                klow = l - 1;
                jchnge = jchnge + 1;
                flag = 0;
            end
        end
        if flag
            l = l - 2;
            while l > klow
                % gee
                c = ad./(cos(2*pi*grid(l))-x);
                err = (c*y'/sum(c) - des(l))*wt(l);
                dtemp = nut*err - comp;
                if dtemp > 0 || jchnge > 0
                    break;
                end
                l = l - 1;
            end
            if l <= klow
                l = iext(j) + 1;
                if jchnge > 0
                    iext(j) = l - 1;
                    j = j + 1;
                    klow = l - 1;
                    jchnge = jchnge + 1;
                else
                    l = l + 1;
                    while l < kup
                        % gee
                        c = ad./(cos(2*pi*grid(l))-x);
                        err = (c*y'/sum(c) - des(l))*wt(l);
                        dtemp = nut*err - comp;
                        if dtemp > 0
                            break;
                        end
                        l = l + 1;
                    end
                    if l < kup && dtemp > 0
                        comp = nut*err;
                        l = l + 1;
                        while l < kup
                            % gee
                            c = ad./(cos(2*pi*grid(l))-x);
                            err = (c*y'/sum(c) - des(l))*wt(l);
                            dtemp = nut*err - comp;
                            if dtemp > 0
                                comp = nut*err;
                                l = l + 1;
                            else
                                break;
                            end
                        end
                        iext(j) = l - 1;
                        j = j + 1;
                        klow = l - 1;
                        jchnge = jchnge + 1;
                    else
                        klow = iext(j);
                        j = j + 1;
                    end
                end
            elseif dtemp > 0
                comp = nut*err;
                l = l - 1;
                while l > klow
                    % gee
                    c = ad./(cos(2*pi*grid(l))-x);
                    err = (c*y'/sum(c) - des(l))*wt(l);
                    dtemp = nut*err - comp;
                    if dtemp > 0
                        comp = nut*err;
                        l = l - 1;
                    else
                        break;
                    end
                end
                klow = iext(j);
                iext(j) = l + 1;
                j = j + 1;
                jchnge = jchnge + 1;
            else
                klow = iext(j);
                j = j + 1;
            end
        end
    end
    while j == nzz
        ynz = comp;
        k1 = min([k1 iext(1)]);
        knz = max([knz iext(nz)]);
        nut1 = nut;
        nut = -nu;
        l = 0;
        kup = k1;
        comp = ynz*1.00001;
        luck = 1;
        flag = 1;
        l = l + 1;
        while l < kup
            % gee
            c = ad./(cos(2*pi*grid(l))-x);
            err = (c*y'/sum(c) - des(l))*wt(l);
            dtemp = err*nut - comp;
            if dtemp > 0
                comp = nut*err;
                j = nzz;
                l = l + 1;
                while l < kup
                    % gee
                    c = ad./(cos(2*pi*grid(l))-x);
                    err = (c*y'/sum(c) - des(l))*wt(l);
                    dtemp = nut*err - comp;
                    if dtemp > 0
                        comp = nut*err;
                        l = l + 1;
                    else
                        break;
                    end
                end
                iext(j) = l - 1;
                j = j + 1;
                jchnge = jchnge + 1;
                flag = 0;
                break;
            end
            l = l + 1;
        end
        if flag
            luck = 6;
            l = ngrid + 1;
            klow = knz;
            nut = -nut1;
            comp = y1*1.00001;
            l = l - 1;
            while l > klow
                % gee
                c = ad./(cos(2*pi*grid(l))-x);
                err = (c*y'/sum(c) - des(l))*wt(l);
                dtemp = err*nut - comp;
                if dtemp > 0
                    j = nzz;
                    comp = nut*err;
                    luck = luck + 10;
                    l = l - 1;
                    while l > klow
                        % gee
                        c = ad./(cos(2*pi*grid(l))-x);
                        err = (c*y'/sum(c) - des(l))*wt(l);
                        dtemp = nut*err - comp;
                        if dtemp > 0
                            comp = nut*err;
                            l = l - 1;
                        else
                            break;
                        end
                    end
                    iext(j) = l + 1;
                    j = j + 1;
                    jchnge = jchnge + 1;
                    flag = 0;
                    break;
                end
                l = l - 1;
            end
            if flag
                flag34 = 0;
                if luck ~= 6
                    iext = [k1 iext(2:nz-nfcns)' iext(nz-nfcns:nz-1)']';
                    jchnge = jchnge + 1;
                end
                break;
            end
        end
    end
    if flag34 && j > nzz,
        if luck > 9
            iext = [iext(2:nfcns+1)' iext(nfcns+1:nz-1)' iext(nzz) iext(nzz)]';
            jchnge = jchnge + 1;
        else
            y1 = max([y1 comp]);
            k1 = iext(nzz);
            l = ngrid + 1;
            klow = knz;
            nut = -nut1;
            comp = y1*1.00001;
            l = l - 1;
            while l > klow
                % gee
                c = ad./(cos(2*pi*grid(l))-x);
                err = (c*y'/sum(c) - des(l))*wt(l);
                dtemp = err*nut - comp;
                if dtemp > 0
                    j = nzz;
                    comp = nut*err;
                    luck = luck + 10;
                    l = l - 1;
                    while l > klow
                        % gee
                        c = ad./(cos(2*pi*grid(l))-x);
                        err = (c*y'/sum(c) - des(l))*wt(l);
                        dtemp = nut*err - comp;
                        if dtemp > 0
                            comp = nut*err;
                            l = l - 1;
                        else
                            break;
                        end
                    end
                    iext(j) = l + 1;
                    jchnge = jchnge + 1;
                    iext = [iext(2:nfcns+1)' iext(nfcns+1:nz-1)' iext(nzz) iext(nzz)]';
                    break;
                end
                l = l - 1;
            end
            if luck ~= 6
                iext = [k1 iext(2:nz-nfcns)' iext(nz-nfcns:nz-1)']';
                jchnge = jchnge + 1;
            end
        end
    end
end

% Inverse Fourier transformation
nm1 = nfcns - 1;
fsh = 1.0e-6;
x(nzz) = -2;
cn = 2*nfcns - 1;
delf = 1/cn;
l = 1;
kkk = 0;
if (edge(1) == 0 && edge(jb) == .5) || nfcns <= 3
    kkk = 1;
end
if kkk ~= 1
    dtemp = cos(2*pi*grid(1));
    dnum = cos(2*pi*grid(ngrid));
    aa = 2/(dtemp - dnum);
    bb = -(dtemp+dnum)/(dtemp - dnum);
end
a = zeros(1,nfcns);
for j = 1:nfcns
    ft = (j-1)*delf;
    xt = cos(2*pi*ft);
    if kkk ~= 1
        xt = (xt-bb)/aa;
        ft = acos(xt)/(2*pi);
    end
    xe = x(l);
    while (xt <= xe) && (xe-xt >= fsh)
        l = l + 1;
        xe = x(l);
    end
    if abs(xt-xe) < fsh
        a(j) = y(l);
    else
        grid(1) = ft;
        % gee
        c = ad./(cos(2*pi*ft)-x(1:nz));
        a(j) = c*y'/sum(c);
    end
    l = max([1, l-1]);
end

dden = 2*pi/cn;

alpha = zeros(1,nfcns);
for j = 1:nfcns
    dnum = (j-1)*dden;
    if nm1 < 1
        alpha(j) = a(1);
    else
        alpha(j) = a(1) + 2*cos(dnum*(1:nm1))*a(2:nfcns)';
    end
end
alpha = [alpha(1) 2*alpha(2:nfcns)]'/cn;
if kkk ~= 1
    p(1) = 2*alpha(nfcns)*bb + alpha(nm1);
    p(2) = 2*aa*alpha(nfcns);
    q(1) = alpha(nfcns-2) - alpha(nfcns);
    for j = 2:nm1
        if j == nm1
            aa = aa/2;
            bb = bb/2;
        end
        p(j+1) = 0;
        sel = 1:j;
        a(sel) = p(sel);
        p(sel) = 2*bb*a(sel);
        p(2) = p(2) + 2*a(1)*aa;
        jm1 = j - 1;
        sel = 1:jm1;
        p(sel) = p(sel) + q(sel) + aa*a(sel+1);
        jp1 = j + 1;
        sel = 3:jp1;
        p(sel) = p(sel) + aa*a(sel-1);
        if j ~= nm1
            sel = 1:j;
            q(sel) = -a(sel);
            q(1) = q(1) + alpha(nfcns - 1 - j);
        end
    end
    alpha(1:nfcns) = p(1:nfcns);
end
if nfcns <= 3
    alpha(nfcns + 1) = 0;
    alpha(nfcns + 2) = 0;
end

alpha=alpha';
% now that's done!

if neg <= 0
    if nodd ~= 0
        h = [.5*alpha(nz-1:-1:nz-nm1) alpha(1)];
    else
        h = .25*[alpha(nfcns) alpha(nz-2:-1:nz-nm1)+alpha(nfcns:-1:nfcns-nm1+2) ...
            2*alpha(1)+alpha(2)];
    end
elseif nodd ~= 0
    h = .25*[alpha(nfcns) alpha(nm1)];
    h = [h .25*[alpha(nz-3:-1:nz-nm1)-alpha(nfcns:-1:nfcns-nm1+3) ...
        2*alpha(1)-alpha(3) ] 0];
else
    h = .25*[alpha(nfcns) alpha(nz-2:-1:nz-nm1)-alpha(nfcns:-1:nfcns-nm1+2) ...
        2*alpha(1)-alpha(2)];
end

%%
function y = remezdd(k, n, m, x)
%REMEZDD Lagrange interpolation coefficients.

%   Author: T. Krauss 1993
%       Was Revision: 1.4, Date: 1994/01/25 17:59:44

y = 1;
q = x(k);
for l=1:m
    xx = 2*(q - x(l:m:n));
    y = y*prod(xx(xx ~= 0));
end
y=1/y;
% EOF


function [b,a]=invfreqz(g,w,varargin)
%INVFREQZ  Discrete filter least squares fit to frequency response data.
%   [B,A] = INVFREQZ(H,W,NB,NA) gives real numerator and denominator
%   coefficients B and A of orders NB and NA respectively, where H is the
%   desired complex frequency response of the system at frequency points W,
%   and W contains the normalized frequency values within the interval [0,
%   Pi] (W is in units of radians/sample).
%
%   INVFREQZ yields a filter with real coefficients. This means that it is
%   sufficient to specify positive frequencies only. The filter fits the
%   data conj(H) at -W, ensuring the proper frequency domain symmetry for a
%   real filter.
%
%   [B,A] = INVFREQZ(H,W,NB,NA,Wt) allows the fit-errors to be weighted
%   versus frequency. The estimates B and A are determined by minimization
%   of sum |B-H*A|^2*Wt over the frequencies in W. LENGTH(Wt) = LENGTH(W) =
%   LENGTH(H).
%
%   [B,A] = INVFREQZ(H,W,NB,NA,Wt,ITER) does another type of fit: Sum
%   |B/A-H|^2*Wt is minimized with respect to the coefficients in B and A
%   by numerical search in at most ITER iterations. The A-polynomial is
%   then constrained to be stable.
%
%   [B,A]=INVFREQZ(H,W,NB,NA,Wt,ITER,TOL) stops the iterations when the
%   norm of the gradient is less than TOL. The default value of TOL is
%   0.01. The default value of Wt is all ones. This default value is also
%   obtained by setting Wt to an empty array [].
%
%   [B,A] = INVFREQZ(H,W,NB,NA,Wt,ITER,TOL,'trace') provides a textual
%   progress report of the iteration.
%
%   [B,A] = INVFREQZ(H,W,'complex',NB,NA,...) creates a complex filter. In
%   this case, no symmetry is enforced and W contains normalized frequency
%   values within the interval [-Pi, Pi].
%
%   % Example:
%   %   Convert a simple transfer function to frequency response data and
%   %   then back to the original filter coefficients. If the system is
%   %   unstable, use the iterative algorithm to find a stable approximation
%   %   to the system.
%
%   b = [1 2 3 2 3];            % Numerator coefficients
%   a = [1 2 3 2 1 4];          % Denominator coefficients
%   [h,w] = freqz(b,a,64);
%   [b1,a1] = invfreqz(h,w,4,5) % a1 has poles outside the unit circle
%   [z,p,k] = tf2zp(b1,a1);     % Get Zero-Pole form
%   fprintf('Stable Approximation to the system:')
%   [b2,a2] = invfreqz(h,w,4,5,[],30) % Stable approximation to system
%   subplot(2,1,1); zplane(b1,a1); title('PZ plot - Unstable system')
%   subplot(2,1,2); zplane(b2,a2); title('PZ plot of stable system')
%
%   See also FREQZ, FREQS, INVFREQS.

%   Author(s): J.O. Smith and J.N. Little, 4-23-86
%              J.N. Little, 4-27-88, revised
%              Lennart Ljung, 9-21-92, rewritten
%              T. Krauss, 10-19-92, trace mode made optional
%   Copyright 1988-2013 The MathWorks, Inc.

% calling sequence is
%function [b,a]=invfreqz(g,w,nb,na,wf,maxiter,tol,pf)
% OR
%function [b,a]=invfreqz(g,w,'complex',nb,na,wf,maxiter,tol,pf)

narginchk(4,9)
if ischar(varargin{1})
    realStr = lower(varargin{1});
    varargin(1) = [];
else
    realStr = 'real';
end
gaussFlag = length(varargin)>3;  % run Gauss-Newton algorithm or not?
if length(varargin)<6
    varargin{6} = [];  % pad varargin with []'s
end
% Checks if 'H' is valid data
signal.internal.sigcheckfloattype(g,'','invfreqz','H');

% Cast to enforce Precision rules
w = signal.internal.sigcasttofloat(w,'double','invfreqz','W','allownumeric');

[nb,na,wf,maxiter,tol,pf] = deal(varargin{:});
% Cast to enforce Precision rules
nb = signal.internal.sigcasttofloat(nb,'double','invfreqz','NB',...
  'allownumeric');
na = signal.internal.sigcasttofloat(na,'double','invfreqz','NA',...
  'allownumeric');
wf = signal.internal.sigcasttofloat(wf,'double','invfreqz','Wt',...
  'allownumeric');
maxiter = signal.internal.sigcasttofloat(maxiter,'double','invfreqz',...
  'ITER','allownumeric');
tol = signal.internal.sigcasttofloat(tol,'double','invfreqz','TOL',...
  'allownumeric');

switch realStr
    case 'real'
        realFlag = 1;
    case 'complex'
        realFlag = 0;
    otherwise
        warning(message('signal:invfreqz:InvalidParam', realStr));
        realFlag = 0;
end

nk=0;T=1; % The code is prepared for arbitrary sampling interval T and for
% constraining the numerator to begin with nk zeros.

nb=nb+nk+1;
if isempty(pf)
    verb=0;
elseif (strcmp(pf,'trace')),
    verb=1;
else
    error(message('signal:invfreqz:NotSupported', pf));
end
if isempty(wf),wf=ones(length(w),1);end
wf=sqrt(wf);

if length(g)~=length(w),error(message('signal:invfreqz:InvalidDimensions', 'H', 'W')),end
if length(wf)~=length(w),error(message('signal:invfreqz:InvalidDimensions', 'Wt', 'W')),end
if any( (w>pi) | (w<0) ) && realFlag
    warning(message('signal:invfreqz:InvalidRegion', 'W', 'INVFREQZ', '''complex'''))
end
[rw,cw]=size(w);    if rw>cw,   w=w';   end
[rg,cg]=size(g);    if cg>rg,   g=g.';  end
[rwf,cwf]=size(wf); if cwf>rwf, wf=wf'; end

nm=max(na,nb+nk-1);
OM=exp(-1i*(0:nm)'*w*T);

%
% Estimation in the least squares case:
%
Dva=(OM(2:na+1,:).').*(g*ones(1,na));
Dvb=-(OM(nk+1:nk+nb,:).');
D=[Dva Dvb].*(wf*ones(1,na+nb));
if realFlag
    R=real(D'*D);
    Vd=real(D'*(-g.*wf));
else
    R=D'*D;
    Vd=D'*(-g.*wf);
end
th=R\Vd;
a=[1 th(1:na).'];b=[zeros(1,nk) th(na+1:na+nb).'];

if ~gaussFlag,return,end

% Now for the iterative minimization
if isempty(maxiter), maxiter = 30; end

if isempty(tol)
    tol = 0.01;
end
indb=1:length(b);indg=1:length(a);
a=polystab(a); % Stabilizing the denominator

% The initial estimate:

GC=((b*OM(indb,:))./(a*OM(indg,:))).';
e=(GC-g).*wf;
Vcap=e'*e; t=[a(2:na+1) b(nk+1:nk+nb)].';
if (verb),
    %messages similar to invfreqs
    clc, disp(['  ' getString(message('signal:invfreqs:INITIALESTIMATE'))]);
    disp([getString(message('signal:invfreqs:CurrentFit')) ' ' num2str(Vcap)])
    disp(getString(message('signal:invfreqs:Parvector')));
    disp(t)
end;

%
% ** the minimization loop **
%
gndir=2*tol+1;l=0;st=0;
while all([norm(gndir)>tol l<maxiter st~=1])
    l=l+1;
    
    %     * compute gradient *
    
    D31=(OM(2:na+1,:).').*(-GC./((a*OM(1:na+1,:)).')*ones(1,na));
    D32=(OM(nk+1:nk+nb,:).')./((a*OM(1:na+1,:)).'*ones(1,nb));
    D3=[D31 D32].*(wf*ones(1,na+nb));
    
    %     * compute Gauss-Newton search direction *
    
    e=(GC-g).*wf;
    if realFlag
        R=real(D3'*D3);
        Vd=real(D3'*e);
    else
        R=D3'*D3;
        Vd=D3'*e;
    end
    gndir=R\Vd;
    
    %     * search along the gndir-direction *
    
    ll=0;k=1;V1=Vcap+1;
    while all([V1 > Vcap ll<20])
        
        t1=t-k*gndir; if ll==19,t1=t;end
        a=polystab([1 t1(1:na).']);
        t1(1:na)=a(2:na+1).';   %Stabilizing denominator
        b=[zeros(1,nk) t1(na+1:na+nb).'];
        GC=((b*OM(indb,:))./(a*OM(indg,:))).';
        V1=((GC-g).*wf)'*((GC-g).*wf); t1=[a(2:na+1) b(nk+1:nk+nb)].';
        if (verb),
            home, disp(int2str(ll))
        end;
        k=k/2;
        ll=ll+1; if ll==20, st=1;end
        if ll==10,gndir=Vd/norm(R)*length(R);k=1;end
    end
    
    if (verb),
        home
        disp(['      ' getString(message('signal:invfreqs:ITERATION')) ' ' int2str(l)])
        disp([getString(message('signal:invfreqs:CurrentFit')) '  ' num2str(V1) '  ' getString(message('signal:invfreqs:PreviousFit')) '  ' num2str(Vcap)])
        disp(getString(message('signal:invfreqs:CurrentParPrevparGNdir')));
        disp([t1 t gndir])
        disp([getString(message('signal:invfreqs:NormOfGNvector')) ' ' num2str(norm(gndir))])
        if st==1,
            disp(getString(message('signal:invfreqs:NoImprovement'))),
            disp(getString(message('signal:invfreqs:IterationsThereforeTerminated'))),
        end
    end
    t=t1; Vcap=V1;
end

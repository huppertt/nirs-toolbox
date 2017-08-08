function [b,a]=invfreqs(g,w,varargin)
%INVFREQS Analog filter least squares fit to frequency response data.
%   [B,A] = INVFREQS(H,W,nb,na) gives real numerator and denominator 
%   coefficients B and A of orders nb and na respectively, where
%   H is the desired complex frequency response of the system at frequency
%   points W, and W contains the frequency values in radians/s.
%   INVFREQS yields a filter with real coefficients.  This means that it is 
%   sufficient to specify positive frequencies only; the filter fits the data 
%   conj(H) at -W, ensuring the proper frequency domain symmetry for a real 
%   filter.
%
%   [B,A]=INVFREQS(H,W,nb,na,Wt) allows the fit-errors to the weighted
%   versus frequency.  LENGTH(Wt)=LENGTH(W)=LENGTH(H).
%   Determined by minimization of sum |B-H*A|^2*Wt over the freqs in W.
%
%   [B,A] = INVFREQS(H,W,nb,na,Wt,ITER) does another type of fit:
%   Sum |B/A-H|^2*Wt is minimized with respect to the coefficients in B and
%   A by numerical search in at most ITER iterations.  The A-polynomial is 
%   then constrained to be stable.  [B,A]=INVFREQS(H,W,nb,na,Wt,ITER,TOL)
%   stops the iterations when the norm of the gradient is less than TOL.
%   The default value of TOL is 0.01.  The default value of Wt is all ones.
%   This default value is also obtained by Wt=[].
%
%   [B,A]=INVFREQS(H,W,nb,na,Wt,ITER,TOL,'trace') provides a textual
%   progress report of the iteration.
%
%   [B,A] = INVFREQS(H,W,'complex',NB,NA,...) creates a complex filter.  In 
%   this case, no symmetry is enforced. 
%
%   % Example:
%   %   Convert a simple transfer function to frequency response data and
%   %   then back to the original filter coefficients. If the system is
%   %   unstable, use invfreqs's iterative algorithm to find a stable 
%   %   approximation to the system.
%
%   b = [1 2 3 2 3];                % Numerator coefficients
%   a = [1 2 3 2 1 4];              % Denominator coefficients
%   [h,w] = freqs(b,a,64);
%   [bb,aa] = invfreqs(h,w,4,5) %  aa has poles in the right half-plane.
%   fprintf('Stable Approximation to the system:')
%   [bbb,aaa] = invfreqs(h,w,4,5,[],30) % stable approximation to system
%
%   See also FREQZ, FREQS, INVFREQZ.

%   Author(s): J.O. Smith and J.N. Little, 4-23-86
%              J.N. Little, 4-27-88, revised
%              Lennart Ljung, 9-21-92, rewritten
%              T. Krauss, 10-22-92, trace mode made optional
%   Copyright 1988-2013 The MathWorks, Inc.
%     

% calling sequence is
%function [b,a]=invfreqs(g,w,nb,na,wf,maxiter,tol,pf)
% OR
%function [b,a]=invfreqs(g,w,'complex',nb,na,wf,maxiter,tol,pf)

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
signal.internal.sigcheckfloattype(g,'','invfreqs','H');
% Cast to enforce Precision rules
w = signal.internal.sigcasttofloat(w,'double','invfreqs','W','allownumeric');

[nb,na,wf,maxiter,tol,pf] = deal(varargin{:});

% Cast to enforce Precision rules
nb = signal.internal.sigcasttofloat(nb,'double','invfreqs','NB',...
  'allownumeric');
na = signal.internal.sigcasttofloat(na,'double','invfreqs','NA',...
  'allownumeric');
wf = signal.internal.sigcasttofloat(wf,'double','invfreqs','Wt',...
  'allownumeric');
maxiter = signal.internal.sigcasttofloat(maxiter,'double','invfreqs',...
  'ITER','allownumeric');
tol = signal.internal.sigcasttofloat(tol,'double','invfreqs','TOL',...
  'allownumeric');

switch realStr
case 'real'
    realFlag = 1;
case 'complex'
    realFlag = 0;
otherwise
    warning(message('signal:invfreqs:InvalidParam', realStr));
    realFlag = 0;
end

nk=0;     % The code is prepared for constraining the numerator to 
	  % begin with nk zeros.

nb=nb+nk+1;
if isempty(pf)
    verb=0;
elseif (strcmp(pf,'trace')),
    verb=1;
else
    error(message('signal:invfreqs:NotSupported', pf));
end
if isempty(wf),wf=ones(length(w),1);end
wf=sqrt(wf);

if length(g)~=length(w)
  error(message('signal:invfreqs:UnmatchedLengths', 'H', 'W'))
end

if length(wf)~=length(w)
  error(message('signal:invfreqs:UnmatchedLengths', 'Wt', 'W'))
end

if any( w(:)<0 ) && realFlag
    warning(message('signal:invfreqs:InvalidWParam', 'W', 'INVFREQS', 'complex'))
end

[rw,cw]=size(w);    if rw>cw,   w=w';   end
[rg,cg]=size(g);    if cg>rg,   g=g.';  end
[rwf,cwf]=size(wf); if cwf>rwf, wf=wf'; end

nm=max(na+1,nb+nk);
indb=nb:-1:1; indg=na+1:-1:1; inda=na:-1:1;

OM=ones(1,length(w));
for kom=1:nm-1
    OM=[OM;(1i*w).^kom]; %#ok<AGROW>
end

%
% Estimation in the least squares case:
%
    Dva=(OM(inda,:).').*(g*ones(1,na));
    Dvb=-(OM(indb,:).');
    D=[Dva Dvb].*(wf*ones(1,na+nb));
    R=D'*D;
    Vd=D'*((-g.*OM(na+1,:).').*wf);
    if realFlag
        R=real(R);
        Vd=real(Vd);
    end
    th=R\Vd;
    a=[1 th(1:na).'];b=[zeros(1,nk) th(na+1:na+nb).'];

if ~gaussFlag,return,end

% Now for the iterative minimization

if isempty(maxiter), maxiter = 30; end
if isempty(tol)
    tol=0.01;
end
% Stabilizing the denominator:
a=apolystab(a,realFlag);

% The initial estimate:

GC=((b*OM(indb,:))./(a*OM(indg,:))).';
e=(GC-g).*wf;
Vcap=e'*e; t=[a(2:na+1) b(nk+1:nk+nb)].';
if (verb),
    % invfreqz using same messages
    clc, disp(['  ' getString(message('signal:invfreqs:INITIALESTIMATE'))]);
    disp([getString(message('signal:invfreqs:CurrentFit')) ' ' num2str(Vcap)])
    disp(getString(message('signal:invfreqs:Parvector')));
    disp(t)
end

%
% ** the minimization loop **
%
gndir=2*tol+1;l=0;st=0;
while all([norm(gndir)>tol l<maxiter st~=1])
    l=l+1;

    %     * compute gradient *

    D31=(OM(inda,:).').*(-GC./((a*OM(indg,:)).')*ones(1,na));
    D32=(OM(indb,:).')./((a*OM(indg,:)).'*ones(1,nb));
    D3=[D31 D32].*(wf*ones(1,na+nb));
    
    %     * compute Gauss-Newton search direction

    e=(GC-g).*wf;
    R=D3'*D3;
    Vd=D3'*e;
    if realFlag
        R=real(R);
        Vd=real(Vd);
    end
    gndir=R\Vd;

    %     * search along the gndir-direction *

    ll=0;k=1;V1=Vcap+1;t1=t;
    while all([V1 > Vcap ll<20]),

        t1=t-k*gndir; if ll==19,t1=t;end
 	    a=[1 t1(1:na).'];
        b=[zeros(1,nk) t1(na+1:na+nb).'];
        a=apolystab(a,realFlag); % Stabilizing the denominator
        t1(1:na)=a(2:na+1).';
        GC=((b*OM(indb,:))./(a*OM(indg,:))).';
        V1=((GC-g).*wf)'*((GC-g).*wf);
        if (verb),
            home, disp(int2str(ll))
        end;
        k=k/2;
        ll=ll+1; if ll==10, gndir=Vd/norm(R)*length(R);k=1;end
        if ll==20,st=1;end
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

function a = apolystab(a,realFlag)
%APOLYSTAB  Stabilize filter, analog
%   inputs: a - denominator polynomial
%           realFlag - 1 for real, 0 for complex
%   returns stabilized denoninator polynomial
if ~isempty(a)
    v=roots(a);
    vind=find(real(v)>0);
    v(vind)=-v(vind);
    a=poly(v);
    if realFlag
        a=real(a);
    end
end

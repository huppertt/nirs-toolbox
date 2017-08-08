function [b,a,b1,b2,sos,g] = maxflat(Z,P,wo,PF)
%MAXFLAT  Maximally flat (a.k.a. generalized Butterworth) digital filter
%design.
%   [B,A] = MAXFLAT(NB,NA,Wn) is a lowpass Butterworth filter with
%   numerator and denominator coefficients B and A of orders NB and NA
%   respectively. Wn is the cut-off frequency at which the filter's
%   magnitude response is equal to 1/sqrt(2) (approx. -3 dB).  Wn must be
%   between 0 and 1, with 1 corresponding to half the sample rate of the
%   filter.  When NA=0, the range of Wn is further restricted due to the
%   smoothness of the filter's response.
%
%   B = MAXFLAT(NB,'sym',Wn) is a symmetric FIR Butterworth filter. NB must
%   be even, and Wn is restricted to a subinterval of [0,1].  The function
%   will raise an error if you specify Wn outside of this subinterval. The
%   filter order, NB, must be less than 140.
%
%   [B,A,B1,B2] = MAXFLAT(NB,NA,Wn) and [B,A,B1,B2] = MAXFLAT(NB,'sym',Wn)
%   return two polynomials B1 and B2 whose product is equal to the
%   numerator polynomial B (i.e., B = CONV(B1,B2)).  B1 contains all the
%   zeros at z=-1, and B2 contains all the other zeros.
%
%   [B,A,B1,B2,SOS,G] = MAXFLAT(NB,NA,Wn) returns the second-order section
%   representation in the matrix SOS and gain G.  For numerical reasons it
%   may be beneficial to use SOS and G in some cases.
%
%   MONITORING THE DESIGN
%      For a textual display of the design table used in the design, use a 
%      trailing 'trace' argument, e.g. MAXFLAT(NB,NA,Wn,'trace').  To
%      display plots of the filter's magnitude, group delay and zeros and
%      poles, use a trailing 'plots' argument, MAXFLAT(NB,NA,Wn,'plots').
%      For both text and plots, use 'both'.
%
%   EXAMPLES:
%
%   % Example 1: IIR design
%                NB = 10; NA = 2; Wn = 0.6; 
%                [b,a,b1,b2] = maxflat(NB,NA,Wn);
%                fvtool(b,a);
%
%   % Example 2: FIR design
%                NB = 10; Wn = 0.6;
%                h = maxflat(NB,'sym',Wn);
%                fvtool(h);
%
%   See also BUTTER, FREQZ, FILTER, DESIGNFILT.

%  by I. W. Selesnick and C. S. Burrus see: Generalized Digital Butterworth
%  Filter Design by I. W. Selesnick and C. S. Burrus, IEEE International
%  Conference on Acoustics, Speech and Signal Processing,
%  May 1996, Volume 3, p. 1367.

%   Copyright 1988-2013 The MathWorks, Inc.

%  required subprograms : spec_table.m, choose.m

TEXT_PF = 0;
PLOT_PF = 0;
SYMM = 0;

% Cast to enforce precision rules.
Z = signal.internal.sigcasttofloat(Z,'double','maxflat','1','allownumeric');
wo = signal.internal.sigcasttofloat(wo,'double','maxflat','3','allownumeric');

if ischar(P)
	SYMM = 1;
	P = 0;  
	if rem(Z,2)==1
	    error(message('signal:maxflat:oddOrder'));
    elseif Z >= 140,
        error(message('signal:maxflat:orderTooLarge'));
	else % (Z is even)
	    Z = Z/2;
	end
end
if nargin > 3  
	if ~ischar(PF)
		error(message('signal:maxflat:invalidPlotFlag'));
	end
	switch lower(PF(1))
	case 't'
		TEXT_PF = 1;
	case 'p'
		PLOT_PF = 1;
	case 'b'
		TEXT_PF = 1;
		PLOT_PF = 1;
	otherwise
		error(message('signal:maxflat:invalidPlotFlag'));
    end
end

P = signal.internal.sigcasttofloat(P,'double','maxflat','2','allownumeric');

SN = 1e-7;
table = spec_table(Z,P,SYMM);
wo = wo*pi;
k = find((table(:,4) < wo/pi+SN) & (table(:,5) > wo/pi-SN), 1, 'last');

%%%%%%%%%%%%%%%%%%%% DISPLAY TABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TEXT_PF
	disp([' ' getString(message('signal:maxflat:Table'))])
	disp(' ')
	disp('    L         M         N         wo_min/pi wo_max/pi') 
	disp(' ')
	disp(table)

	if SYMM
	disp(getString(message('signal:maxflat:LAndMAreDoubledForSymmetricFIRFilters')))
	disp(' ')
	end

	%  L+M : total number of zeros
	%  N   : total number of poles
	%
	%  L : number of zeros at z=-1
	%  M : number of zeros contributing to passband flatness
	%  N : number of poles
end

if PLOT_PF
        clf
        set(gcf,'nextplot','add');
        a1 = axes('units','normalized','position',[.1 .6 .8 .3]);
        a2 = axes('units','normalized','position',[.1 .1 .35 .4]);
        a3 = axes('units','normalized','position',[.55 .1 .35 .4]);
end

%%%%%%%%%%%%%%%%%%%% error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(k)
    error(message('signal:maxflat:invalidCutoff', sprintf( '%g', (1 + SYMM)*Z ), sprintf( '%g', P ), sprintf( '%g', table( 1, 4 ) ), sprintf( '%g', table( end, 5 ) )));
end

%%%%%%%%%%%%%%%%%%%% CALCULATE x-Domain function %%%%%%%%%%%%%%%%%%
L = table(k,1);
M = table(k,2);
N = table(k,3);
xo = (1-cos(wo))/2;
%Fo = (1/2)^2;
Fo = 1/2;
if SYMM
    Fo = sqrt(Fo);
end

if N == 0
	%%%%%%%%%%%% FIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	k = M-1:-1:0;
	R = [0, choose(M-k-1,0).*choose(L+k-1,k)];
   	T = [choose(M-k-2,-1).*choose(L+k,k), 0];
        c1 = (Fo/(1-xo)^L - polyval(R,xo));
        c2 = polyval(T,xo);
        S = c2*R + c1*T;
	Q = 1;
elseif M == 0
	%%%%%%%%%%%% ALL ZEROS at z = -1 %%%%%%%%%%%%%%%%%%%%%%%%%%
	k = min([L,N]):-1:0;
	Q = choose(L,k).*(-1).^k;
	if L < N
		Q = [zeros(1,N-L), Q];
	end
	c1 = ((1/Fo)*(1-xo)^L - polyval(Q,xo))/(xo^N);
	Q(1) = Q(1) + c1;
	S = 1;
else % M > 0
	%%%%%%%%%%%% Some Zeros contribute to passband %%%%%%%%%%%%
	k = M-1:-1:0;
	R = [0, choose(M+N-k-1,N).*choose(L-N+k-1,k)];
	T = [choose(M+N-k-2,N-1).*choose(L-N+k,k), 0];
	k = N:-1:0;
	AGR = choose(M+N-k-1,M-1).*choose(M+L-1,k).*(-1).^k;
	AGT = choose(M+N-k-1,M-1).*choose(M+L-1,k-1).*(-1).^(k-1);
        c1 = (((1-xo)^L)*polyval(R,xo) - Fo*polyval(AGR,xo));
        c2 = (Fo*polyval(AGT,xo) - ((1-xo)^L)*polyval(T,xo));
        S = c2*R + c1*T;
        Q = c2*AGR + c1*AGT;
end

%%%%%%%%%%%%%%%%%%%% convert to z-domain %%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = 1 - 2*roots(S);
br1 = tmp+sqrt(tmp.^2-1);
br2 = tmp-sqrt(tmp.^2-1);
br = sort([br1; br2]+eps*sqrt(-1));     % sort by absolute value
if ~SYMM
    if c1 ~= 0
        br = br(1:M);                   % take roots INSIDE unit circle
    else
        br = [br(1:M-1); 0];            % take roots INSIDE unit circle
    end
end
b2 = real(poly(br));

Qrts = roots(Q);
ar1 = 1-2*Qrts+sqrt((1-2*Qrts).^2-1);
ar2 = 1-2*Qrts-sqrt((1-2*Qrts).^2-1);
ar = sort([ar1; ar2]+eps*sqrt(-1));     % sort by absolute value
if ~SYMM
    ar = ar(1:N);                       % take roots INSIDE unit circle
end
a = real(poly(ar));

%%%%%%%%%%%%%%%%%%%% construct b1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SYMM
   L = 2*L;
end

b1 = 1;
b = b2;
for t = 1:L
   b1 = conv(b1,[1 1]);
   b = conv(b,[1 1]);
end

%%%%%%%%%%%%%%%%%%% build sos matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sos,g] = buildsos(br,L,ar,a,b2);
%%%%%%%%%%%%%%%%%%%% normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(sum(b2)) < 100*eps
   % terminate program
   error(message('signal:maxflat:SignalErr'))
else
   b2 = b2*sum(a)/(sum(b1)*sum(b2));
   b = conv(b1,b2);
end

%%%%%%%%%%%%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(a) && PLOT_PF

    % ---- PLOT MAGNITUDE ----
    [Mb1] = h2mag(b1);
    [Mb2] = h2mag(b2);
    [Ma,w]  = h2mag(a);
    Mag = Mb1.*Mb2./Ma;
    axes(a1)
    plot(w/pi,Mag)
    xlabel('\omega/\pi')
    ylabel(getString(message('signal:maxflat:Magnitude')))
    title(getString(message('signal:maxflat:FrequencyResponse')))
    axis([0 1 -.1 1.1])

    % ---- PLOT POLE-ZERO DIAGRAM ----
    I = sqrt(-1);
    SN = 1e-7;
    axes(a2)
    circ = exp((0:100)/100 * 2 * pi * I);
    plot(circ,'--')
    hold on
    if length(ar)
        plot(ar+SN*I,'x')
    end
    if length(br)
        plot(br+SN*I,'o')
    end
    plot(-1+SN*I,'o')
    hold off
    axis([-1 1 -1 1]*1.2)
    xlabel(getString(message('signal:maxflat:Real')))
    ylabel(getString(message('signal:maxflat:Imaginary')))
    if SYMM
       title(getString(message('signal:maxflat:ZeroPlotzerosOutsideCircleNotShown')))
    else
       title(getString(message('signal:maxflat:PolezeroPlot')))
    end
    s1 = ['<- deg ',num2str(L)];
    text(-0.9,0,s1);
    axis('square')

    % ---- PLOT GROUP DELAY ----
    gd = grpdelay(b2,a,w)+L/2;
    axes(a3)
    plot(w/pi,gd)
    mingd = min(gd);
    maxgd = max(gd);
    if (maxgd-mingd)<(1e-10)*(abs(maxgd)+abs(mingd))
        set(a3,'ylim',[0 2*maxgd])
    end
    xlabel('\omega/\pi')
    ylabel(getString(message('signal:maxflat:Samples')))
    title(getString(message('signal:maxflat:GroupDelay')))
    set(gcf,'nextplot','replace');

end


function table = spec_table(Z,P,SYMM)
%
% table = spec_table(Z,P,SYMM);
%
% General Butterworth Filter Design
% This program decides how many zeros should lie at 
% z = -1 and how many should contribute to the passband
%
% Z : total number of zeros
% P : total number of poles
%
% % Example:
%
%   Z = 10; P = 2; 
%   table = spec_table(Z,P);
%
%  Ivan Selesnick, Rice University,
%  September 1995.
%  see: Generalized Digital Butterworth Filter Design

%  required subprograms : choose.m

if (Z <= P) 
	table = [Z 0 P 0 1];
	return
else
	table = [];
end

%Fo = (1/2)^2;   % value of mag squared at half-mag freq.
Fo = 1/2;   % value of mag squared at half-mag freq.
if SYMM
    Fo = sqrt(Fo);
end

%%%%%%% begin with all zeros at z = -1 %%%%%%%
L = Z;
N = P;
if N > 0
	if rem(N,2) == 0
		c = 0;			% N even
	else
		c = choose(L-1,N);	% N odd
	end
	%%%%%%% construct polynomial for checking boundary %%%%%%
	Y = zeros(1,L+1);
	for i = 0:N
		Y(i+1) = choose(L,i)*(1-Fo)*(-1)^i;
	end
	for i = N+1:L
		Y(i+1) = choose(L,i)*(-1)^i;
	end
	Y(N+1) = Y(N+1) - c*Fo;
	Y = Y(L+1:-1:1);
	%%% extract appropriate root %%%
	r = specroot(Y);
	wo_max = acos(1-2*r);

	table = [Z, 0, P, 0, wo_max/pi];
end

%%%%%%%% first order case %%%%%%
if (Z==1) && (P==0)
    % In this case, we cannot use second generalization because it requires
    % Z > P+1, so we can only use first generalization. Therefore, using
    % Selesnick's notation, L = 1, N = 0, M = 0, F(x) = (1-x)
    % In this section we choose boudary value c=0. If we first calcuate c
    % using the cutoff frequency and then try to come up the limit for wo,
    % we will have the same limit.
    
    r = 1-Fo;
    wo_max = acos(1-2*r);
    
    table = [1, 0, 0, 0, wo_max/pi];
end

%%%%%%%%% go through each and check %%%%%%%%%
for M = 1:Z-P-1
	L = Z-M;
%     k = M-1:-1:0;
%     R = [0, choose(M+N-k-1,N).*choose(L-N+k-1,k)];
%     T = [choose(M+N-k-2,N-1).*choose(L-N+k,k), 0];
	k = M+L:-1:0;
        GR = choose(M+N-k-1,M-1).*choose(M+L-1,k).*(-1).^k;
        GT = choose(M+N-k-1,M-1).*choose(M+L-1,k-1).*(-1).^(k-1);
	AGR = GR(M+L+1-N:M+L+1);
	AGT = GT(M+L+1-N:M+L+1);

	%%% ------- wo_max ---------------------------------- %%%
        if rem(N,2) == 0
                c = (L-N)/(M+N);  % N even
        else
                c = (L-N)/N;      % N odd
        end
        Y = [zeros(1,M+L-N), Fo*AGR+c*Fo*AGT]-GR-c*GT;
	%%%%%%% extract appropriate root %%%
        r = specroot(Y);
        if isempty(r)
          error(message('signal:maxflat:UnableToDesign'))
        end        
        wo_max = acos(1-2*r);

        %%% ------- wo_min ---------------------------------- %%%
        if rem(N,2) == 0
                c = -1;		% N even
		Y = [zeros(1,M+L-N), Fo*AGR+c*Fo*AGT]-GR-c*GT;
        else
%             c = inf; % N odd
		Y = [zeros(1,M+L-N), Fo*AGT]-GT;
        end
        %%%%%%% extract appropriate root %%%
        r = specroot(Y);
        if isempty(r)
          error(message('signal:maxflat:UnableToDesign'))
        end
        
        wo_min = acos(1-2*r);
        
	table = [table; [L,M,N, wo_min/pi, wo_max/pi]];
end

if N > 0
	M = Z-P;
	L = Z-M;
	table = [table; [L,M,N, table(Z-P,5), 1]];
end

function a = choose(n,k)
%
% a = choose(n,k)
% BINOMIAL COEFFICIENTS
%
% allowable inputs:
%   n : integer, k : integer
%   n : integer vector, k : integer
%   n : integer, k : integer vector
%   n : integer vector, k : integer vector (of equal dimension)
%

nv = n;
kv = k;
if (length(nv) == 1) && (length(kv) > 1)
	nv = nv * ones(size(kv));
elseif (length(nv) > 1) && (length(kv) == 1)
	kv = kv * ones(size(nv));
end
a = nv;
for i = 1:length(nv)
   n = nv(i);
   k = kv(i);
   if n >= 0
      if k >= 0
         if n >= k
            c = prod(1:n)/(prod(1:k)*prod(1:n-k));
         else
            c = 0;
        end
     else
        c = 0;
     end
   else
      if k >= 0
         c = (-1)^k * prod(1:k-n-1)/(prod(1:k)*prod(1:-n-1));
      else
         if n >= k
            c = (-1)^(n-k)*prod(1:-k-1)/(prod(1:n-k)*prod(1:-n-1));
         else
            c = 0;
         end
      end
   end
   a(i) = c;
end

function r = specroot(Y)
r = roots(Y);
r = r(imag(r)==0);
[temp,i] = min(abs(r-0.5)); %#ok
r = r(i);
%  other alternatives:
%SN = 1e-5;
%r = newton(Y,[SN 1-SN]);
%f = @(x,p) polyval(p,x);
%r = fzero(f,.5,optimset('TolX',eps,'Display','off'),Y);


function [M,w] = h2mag(h)
% [M,w] = h2mag(h)
L = 2^9;
f = fft(h,2*L);
M = abs(f);
M = M(1:L+1);
w = (0:L)'*pi/L;

%------------------------------------------------------------------
function [sos,g] = buildsos(br,L,ar,a,b2)

% Construct numerator portion, only use first 3 columns
if ~isempty(br),
    tempzeroimagidx = find(abs(imag(br)) <= 2^-32);   %detect zero imaginary part
    br(tempzeroimagidx) = real(br(tempzeroimagidx));  %force real number to be real
    [sosnum,gnum] = zp2sos(br,zeros(size(br)),1); %#ok
else
    sosnum = [];
    gnum = 1; %#ok
end
% Add the zeros at -1
sosnumatm1 = repmat([1/4 1/2 1/4 0 0 0],floor(L/2),1);
if rem(L,2),
    % Add one more zero
    sosnumatm1(end+1,:) = [1/2 1/2 0 0 0 0];
end
sosnum = [sosnum;sosnumatm1];

if ~isempty(ar),
    % Remove small imaginary parts
    indx = find(abs(imag(ar)) < eps^(2/3) & ...
        (real(ar) == 0 | abs(imag(ar)./real(ar)) < sqrt(eps))); 
    ar(indx) = real(ar(indx));
    [sosden,gden] = zp2sos(zeros(size(ar)),ar,1); %#ok
else
    sosden = repmat([1 0 0 1 0 0],size(sosnum,1),1);
    gden = 1; %#ok
end

% Make sure sizes are the same before concatenating
sosnum = [sosnum;repmat([1 0 0 1 0 0],size(sosden,1)-size(sosnum,1),1)];
sosden = [sosden;repmat([1 0 0 1 0 0],size(sosnum,1)-size(sosden,1),1)];

sos = [sosnum(:,1:3),sosden(:,4:6)];
g = sum(a)/sum(b2);

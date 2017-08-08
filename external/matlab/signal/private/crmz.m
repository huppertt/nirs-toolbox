function [h,a,delta,not_optimal,iext,HH,EE, ...
          M_str,HH_str,h_str,Lf,Lb,Ls,Lc,A] = ...
          crmz( L, sym, Lfft, indx_edges, PLOTS, TRACE )
%CRMZ Complex Remez multiple-exchange filter design algorithm.
%   Designs FIR filters with arbitrary magnitude and 
%   phase specifications; reduces to the Remez (Parks-McClellan)
%   algorithm in the real case.

%   Authors: L. Karam, J. McClellan
%   Revised: 22-Oct-96, D. Orofino
%
%   Copyright 1988-2004 The MathWorks, Inc.

global DES_CRMZ WT_CRMZ GRID_CRMZ TGRID_CRMZ IFGRD_CRMZ

J = 1i;
N1 = 0;

if nargin==1 & isstruct(L),
 if isstr(L.plot) & strcmp(L.plot,'plot-result'),
    N2 = N1 + L.L - 1;
    H = L.HH .* exp(-J*2*pi*TGRID_CRMZ*(N1+N2)/2);
    crmz_plresult( H, L.EE, L.iext, L.iext, L.indx_edges, N1, N2, L.delta, L.sym );
    drawnow
 end
 return
end

if TRACE,
  disp(' ');
  disp('             ***********************************************');
  disp(['             ' getString(message('signal:crmz:ComplexRemez'))]);
  disp('             ***********************************************');
end

if nargin<2, sym = 0; end;

N2 = N1 + L - 1;   % Last index of impulse response 
is_odd = rem(L,2); 
if (is_odd)
  Lf = (L+1)/2;
  Ws = ['W(:,2:Lf)'];
  hr_str = ['[a([Lc:-1:2])/2; a(1); a([2:Lc])/2]'];
  hi_str = ['[-a([Lb:-1:(Lc+1)])/2;    0; a([(Lc+1):Lb])/2]'];
  ph_str = ['ones(length(TGRID_CRMZ),1)'];
else
  Lf = L / 2;
  Ws = ['W'];
  hr_str = ['[a([Lc:-1:1])/2; a([1:Lc])/2]'];
  hi_str = ['[-a([Lb:-1:(Lc+1)])/2; a([(Lc+1):Lb])/2]'];
  ph_str = ['exp(-J*pi*TGRID_CRMZ)'];
end;
Lc = Lf;
Ls = L - Lf;

switch sym
case 0,
  Lb = L;
  M_str = ['[cos(W)  sin(' Ws ')]'];
  h_str = [hr_str '+J*' hi_str];
  HH_str = ['HH = fftshift( fft(hc,Lfft) ).*' ph_str ';'];
case 1,
  % Even-symmetry:
  Lb = Lc;
  M_str = ['cos(W)'];
  h_str = [hr_str];
  HH_str = ['HH=fft(hc,Lfft); HH((Lfft/2+2):Lfft)=[]; HH=HH.* ' ph_str ';'];
case 2,
  % Odd-symmetry:
  Lb = Ls;
  Lc = 0;
  M_str = ['sin(' Ws ')'];
  h_str = ['J*' hi_str];
  HH_str = ['HH=fft(hc,Lfft); HH((Lfft/2+2):Lfft)=[]; HH=HH.* ' ph_str ';'];
end

%------------------------
A = DES_CRMZ .* exp(J*2*pi*GRID_CRMZ*(N1+N2)/2);
if PLOTS
  % clf;
  plot(GRID_CRMZ*2,abs(DES_CRMZ),'--')
  hold on
  plot(GRID_CRMZ*2, angle(DES_CRMZ),'-')
  hold off 
  title(getString(message('signal:crmz:DesiredFilterFrequencyResponse')))
  xlabel(getString(message('signal:crmz:NormalizedFrequency'))) 
  ylabel(getString(message('signal:crmz:MagnitudedashedAndPhasesolid')))
  drawnow
end

%-----------------------
if TRACE, disp(getString(message('signal:crmz:CalculatingInitialSolution'))), end

vec_edges    = TGRID_CRMZ(indx_edges);
edges        = reshape(vec_edges,2,length(vec_edges)/2)';
[fext, iext] = crmz_guess(edges, GRID_CRMZ, Lb);

it = 0; delta = 0.0; delta_old = -1; no_stp = 1;

exactTol = eps^(2/3);  % tolerance for exactness (arbitrary)
% if the maximum error relative to the absolute maximum of the
%  desired function is less than exactTol, the filter designed
%  is considered 'exact' and the iteration is terminated successfully.
last_ee=zeros(1,10);
while (no_stp)
  it = it + 1;
  delta_old = abs(delta);
  fext = GRID_CRMZ(iext);
  W    = 2*pi*fext*([0:(Lf-1)]+(~is_odd)*0.5);
  Mb   = eval(M_str);
  M    = [ Mb ((-1).^[0:Lb]')./WT_CRMZ(iext) ];
  a    = M \ A(iext);
  delta = a(Lb+1);
  h  = eval(h_str);
  hc = crmz_rotate(zeropad(h,Lfft-L), -Ls);
  eval(HH_str);
  W  = 2*pi*vec_edges*([0:(Lf-1)]+(~is_odd)*0.5);
  Mb = eval(M_str);

  HH(indx_edges) = Mb * a(1:Lb);
  EE = WT_CRMZ .* (A - HH(IFGRD_CRMZ));
  EE(iext)   = delta*( (-1) .^ [2:length(iext)+1] ); 
  [jext,EEj] = crmz_find(EE,iext);  
  if TRACE,
    fprintf('(%3.0f): ', it) 
    fprintf('delta = %8.5f+j%8.5f, ', real(delta), imag(delta)) 
    fprintf('|delta| = %5.2e, ', abs(delta));
    fprintf('e_max = %5.2e\n', max(abs(EE)));
  end
  if PLOTS,
    crmz_plerror(EE,HH,iext,jext,indx_edges,GRID_CRMZ,delta,sym);
    drawnow
  end
  e_max=max(abs(EE));
  last_ee=[e_max last_ee(1:end-1)];
  s=max(abs(last_ee-e_max(:,ones(1,10))));
  if all(iext == jext') | (max(abs(EE))/max(abs(A))<exactTol)...
          | ((s<exactTol) & (it>10)),
    no_stp = 0;
  end
  iext = jext';
end
%%% end of Complex Remez algorithm %%%%

%%%%% Assessing Optimality of the Solution %%%%
tlr = abs(delta)/100;  %% Needed due to limited machine accuracy %% 
if (e_max <= (abs(delta)+tlr)) | (e_max/max(abs(A))<exactTol)
  not_optimal = 0; 
  if TRACE,
    fprintf(getString(message('signal:crmz:OptimalSolutionObtained')))
    fprintf('(%3.0f): ', it) 
    fprintf('delta = %8.5f+j%8.5f, ', real(delta), imag(delta)) 
    fprintf('|delta| = %5.2e, ', abs(delta));
    fprintf('e_max = %5.2e\n',e_max);
  end
else
  not_optimal = 1; 
  if TRACE
    crmz_chckopt( EE, delta, tlr, Lfft )
    fprintf('(%3.0f): ', it) 
    fprintf('delta = %8.5f+j%8.5f, ', real(delta), imag(delta)) 
    fprintf('|delta| = %5.2e, ', abs(delta));
    fprintf('e_max = %5.2e\n',e_max);
  end
end

%%%%%%%%%%%%%  End Complex Remez Stage  %%%%%%%%%%%%%%%%%%%%%%%%

%========================================================

function crmz_chckopt( EE, delta, tlr, Lfft )
% Checks degree of suboptimality of solution not optimal 
% on desired frequency bands. 
%
subint = find( abs(EE) <= (abs(delta)+tlr) );
rel_size = round((length(subint)/Lfft)*100);
disp(getString(message('signal:crmz:SolutionInfo',sprintf('%3.0f',rel_size))));

%========================================================

function  [fnew,enew] = crmz_find( error, fold )
%  Construct new subset of points using exchange rules

%    Authors: LJ Karam and JH McClellan
%  Reference: "Complex Chebyshev approximation for FIR filter design",
%             IEEE Trans. on Circuits and Systems II, March 1995. 

fold      = fold(:);
Nx        = length( fold );
Ngrid     = length( error );
if error(fold(1))~=0
    sgn_error = error(fold(1))/abs(error(fold(1)));
else
    sgn_error = 1;
end
error     = real( conj(sgn_error)*error );
delta     = min( abs( error(fold) ) ); %--- present value of delta
up        = sign(error(fold(1))) > 0;
if up,
  fence = [ [1;fold(1:2:Nx)]  [fold(1:2:Nx);Ngrid] ];
else
  fence = [ [1;fold(2:2:Nx)]  [fold(2:2:Nx);Ngrid] ];
end
Lf  = length(fence(:,1));
emn = zeros(Lf,1); imn = emn; emx = emn; imx = emn;
for i=1:Lf,
   [emn(i),imn(i)] = min( error(fence(i,1):fence(i,2) ));
   imn(i) = imn(i) + fence(i,1) -1;
end
imn(find((emn > -delta))) = [];
fence = [ [1;imn]  [imn;Ngrid] ];
Lf = length(fence(:,1));
for i=1:Lf
   [emx(i),imx(i)] = max( error(fence(i,1):fence(i,2) ));
   imx(i) = imx(i) + fence(i,1) -1;
end
imx(find((emx < delta))) = [];
fnew = sort([imx;imn]);
Nf   = length(fnew);
if ( Nf > Nx)
 if ( abs(error(fnew(1))) >= abs(error(fnew(Nf))) )
   fnew = fnew(1:Nx);
 else
   fnew = fnew(Nf-Nx+1:Nf);
 end
end
fnew = fnew';
enew = error(fnew);

%========================================================

function  [fext,iext] = crmz_guess( edges, grid, nfcns )
%CRMZ_GUESS   generate initial guess of "EXTREMAL FREQS"
%               for firpm exchange algorithm
%   usage:
%        [Fext,Iext] = crmz_guess( Edges, Grid, Nfcns )
%
%     Edges:  band-edges moved onto the grid (Revised edges)
%     Grid:   frequencies on the fine grid
%     Nfcns:  number of basis functions in the approx
%     Fext:   initial guess of "extremal frequencies"
%     Iext:   indices for the "extremal frequencies"
%
TOL    = 5*eps;
next   = nfcns+1;
Nbands = length(edges(:,2));
tt     = edges;
merged = tt;
if Nbands > 1,
  jkl = find( abs(tt(1:(Nbands-1),2)-tt(2:Nbands,1)) > TOL );
  merged = [ [tt(1,1); tt(jkl+1,1)] , [tt(jkl,2); tt(Nbands,2)] ];
end
Nbands = length(merged(:,1));
bw = merged(:,2)-merged(:,1);
if any(bw < 0),
   edges
   error(message('signal:crmz:InternalErrorNegBW'));
end
percent_bw = bw / sum(bw);
fext = zeros(next,1);
n = 0;  i = 1;
while (n < next) & (i <= Nbands),
   nfreqs_i = min( next-n, ceil( percent_bw(i)*next ) );
   if( nfreqs_i == 0 )
      n=n+1;   fext(n) = merged(i,1);
   else
      fext(n+[1:nfreqs_i]) = ...
           linspace(merged(i,1),merged(i,2),nfreqs_i);
      n = n+nfreqs_i;
   end
   i = i+1;
end
iext = 0*fext;
for i=1:next
   [tt,iext(i)] = min( abs(fext(i)-grid) );
end
if any(diff(iext) == 0),
   error(message('signal:crmz:InternalErrorGridPoint'))
end
fext = grid(iext);

%========================================================

function rotated = crmz_rotate(x,num_places)
%CRMZ_ROTATE     circular shift of columns in a matrix
%    crmz_rotate(V,r)
%        circularly shifts the elements in the columns of V
%        by r places right (r>0); or r places left (r<0).
%           (Right is down; left is up.)
%        If the input is a row or column vector, the shift is
%        performed on the vector.
%        If the input is a signal matrix, each column is shifted
%
[M,N] = size(x);
if M > 1,
   % rotate columns
   num_places = mod(num_places,M);  % make num_places in range [0,M-1]
   rotated    = [ x(M-num_places+1:M,:); x(1:M-num_places,:) ];
else  % Assume N > 1
   % rotate row vector
   num_places = mod(num_places,N);  % make num_places in range [0,N-1]
   rotated    = [ x(N-num_places+1:N) x(1:N-num_places) ];
end

%========================================================

function crmz_plerror2( H, EE, iext, jext, indx_edges, N1, N2, delta, sym )
% crmz_plerror2
% plot error magnitude and real phase-rotated error. 

global DES_CRMZ WT_CRMZ GRID_CRMZ TGRID_CRMZ IFGRD_CRMZ

J = 1i;

EE1 = EE(iext(1));
EE_pl = [];
grd_pl = [];
i = 1; l = 0;
while (i < (length(indx_edges)-1))
  lo = l + 1;
  l  = l + length([indx_edges(i):indx_edges(i+1)]);
  grd_pl = [grd_pl;GRID_CRMZ(lo:l);nan];
  EE_pl = [EE_pl;EE(lo:l);nan];
  i = i + 2;
end;
lo = l + 1;
l  = l + length([indx_edges(i):indx_edges(i+1)]);
grd_pl = [grd_pl;GRID_CRMZ(lo:l)];
EE_pl  = [EE_pl;EE(lo:l)];
EEE_pl = [real(EE_pl * conj(EE1/abs(EE1)) )-3*abs(delta),abs(EE_pl)]; 
EEE    = [real(EE * conj(EE1/abs(EE1)) )-3*abs(delta), abs(EE)]; 
delr   = real(delta);   deli = imag(delta);
axis('auto')
plot( grd_pl*2, EEE_pl, '-', [GRID_CRMZ(1), 0.5]*2, abs(delta)*[1;1], '--',...
 [GRID_CRMZ(1), 0.5]*2, abs(delta)*[1 -1;1 -1]-3*abs(delta), ':',...
  GRID_CRMZ(jext)*2, EEE(jext,:), 'o', GRID_CRMZ(iext)*2, EEE(iext,:), '+')
V = axis;
axis([-1*(~sym) 1 V(3) V(4)])

%========================================================

function crmz_plresult( H, EE, iext, jext, indx_edges, N1, N2, delta, sym )
% crmz_plresult
%
% Generates subplots of final result with
%    1. Magnitude response
%    2. Passband Phase 
%    3. Error magnitude and Real rotated error
%    4. Phase error in passband 
% Uses grp_delay.m

global DES_CRMZ WT_CRMZ GRID_CRMZ TGRID_CRMZ IFGRD_CRMZ

J = 1i;

% Compute EE_pl for plotting error in bands:
EE_pl = []; grd_pl = []; i = 1; l = 0;
while (i < (length(indx_edges)-1))
  lo = l + 1;
  l = l + length([indx_edges(i):indx_edges(i+1)]);
  grd_pl = [grd_pl;GRID_CRMZ(lo:l);nan];
  EE_pl = [EE_pl;EE(lo:l);nan];
  i = i + 2;
end
lo = l + 1;
l  = l + length([indx_edges(i):indx_edges(i+1)]);
grd_pl = [grd_pl;GRID_CRMZ(lo:l)];
EE_pl = [EE_pl;EE(lo:l)];
E_pl = EE_pl .* exp(-2i*pi*grd_pl*(N1+N2)/2);  % Needed for plot 

hfig = figure;
subplot(111)
axis('auto');
subplot(221)
dbMax = 20*log(max(abs(H(:)))); 
%dbMin = 20*log(min(abs(H(:))));
%dbRange = dbMax - dbMin; 
dbRange = 80;
plot(TGRID_CRMZ*2, db(abs(H(:)),dbRange,dbMax), '-') 
ph = gca;
set(ph,'box','off');
xlabel(getString(message('signal:crmz:NormalizedFrequency'))); ylabel(getString(message('signal:crmz:MagnitudelogScale')))
subplot(222);
if length(indx_edges) > 2
  pas_edges = TGRID_CRMZ( indx_edges(3:4) );
else
  pas_edges  = TGRID_CRMZ( indx_edges(1:2) );
end
pas_edges = pas_edges';
pas_edges = pas_edges(:);
i = 1; 
grd_pl = [];
aH_pl = [];
aD_pl = []; 
while (i < (length(pas_edges)-1))
   psb = find((GRID_CRMZ >= pas_edges(i)) & (GRID_CRMZ <= pas_edges(i+1)));
   grd_pl = [grd_pl;GRID_CRMZ(psb);nan];
   aH_pl = [aH_pl;angle(H(IFGRD_CRMZ((psb))));nan];
   aD_pl = [aD_pl;angle(DES_CRMZ(psb));nan];
   i = i + 2; 
end
psb = find((GRID_CRMZ >= pas_edges(i)) & (GRID_CRMZ <= pas_edges(i+1)));
grd_pl = [grd_pl;GRID_CRMZ(psb)];
aH_pl = [aH_pl;angle(H(IFGRD_CRMZ((psb))))];
aD_pl = [aD_pl;angle(DES_CRMZ(psb))];
plot(grd_pl*2, aH_pl, '-')
hold on
plot(grd_pl*2, aD_pl, '--')
ph = gca;
set(ph,'box','off');
xlabel(getString(message('signal:crmz:NormalizedFrequency'))), ylabel(getString(message('signal:crmz:PassbandPhaseradians')))
V = axis;
V(1) = pas_edges(1);
V(2) = pas_edges(length(pas_edges));
axis([2*V(1) 2*V(2) V(3) V(4)]);
x = (V(1)+V(2))/(2.0);
text(2*x,V(4),'( -- ideal )')
hold off
subplot(224)
e_angle = mod(aD_pl,2*pi)-mod(aH_pl,2*pi);
plot(grd_pl*2,e_angle)
xlabel(getString(message('signal:crmz:NormalizedFrequency'))), ylabel(getString(message('signal:crmz:PassbandPhaseError')))
ph = gca;
set(ph,'box','off');
V = axis;
V(1) = pas_edges(1); 
V(2) = pas_edges(length(pas_edges)); 
axis([2*V(1) 2*V(2) V(3) V(4)]);
subplot(223)
crmz_plerror2( H, EE, iext, jext, indx_edges, N1, N2, delta, sym )
ph = gca;
set(ph,'box','off');
axis('off')
text(.95,-4.5*abs(delta),'1')
text(-0.0108,-4.5*abs(delta),'0')
if (~sym)
 text(-0.52*2,-4.5*abs(delta),'-1')
 str = sprintf('%5.4f',abs(delta));
 text(-0.75*2,abs(delta),str,'fontsize',10); 
 text(-0.75*2,-2*abs(delta),str,'fontsize',10);
 str = sprintf('%5.4f',-abs(delta));
 text(-0.78*2,-4*abs(delta),str,'fontsize',10);
 text(-0.22*2,1.6*abs(delta),getString(message('signal:crmz:ErrorMagnitude')),'fontsize',10);
 text(-0.25*2,-1.5*abs(delta),getString(message('signal:crmz:RealRotatedError')),'fontsize',10);
else 
 str = sprintf('%5.4f',abs(delta));
 text(-0.14*2,abs(delta),str,'fontsize',10); 
 text(-0.14*2,-2*abs(delta),str,'fontsize',10);
 str = sprintf('%5.4f',-abs(delta));
 text(-0.16*2,-4*abs(delta),str,'fontsize',10);
 text(0.155*2,1.6*abs(delta),getString(message('signal:crmz:ErrorMagnitude')),'fontsize',10);
 text(0.125*2,-1.5*abs(delta),getString(message('signal:crmz:RealRotatedError')),'fontsize',10);
end 
xlabel(getString(message('signal:crmz:NormalizedFrequency')))

%========================================================

function crmz_plerror(EE,HH,iext,jext,indx_edges,grd,delta,sym)
%crmz_plerror
% plot error magnitude, real phase rotated error, 
%      error real part, error imaginary part. 
%
clf;
EE1 = EE(iext(1));
EE_pl = []; grd_pl = []; i = 1; l = 0;
while (i < (length(indx_edges)-1))
  lo = l + 1;
  l  = l + length([indx_edges(i):indx_edges(i+1)]);
  grd_pl = [grd_pl;grd(lo:l);nan];
  EE_pl  = [EE_pl;EE(lo:l);nan];
  i = i + 2;
end;
lo = l + 1;
l  = l + length([indx_edges(i):indx_edges(i+1)]);
grd_pl = [grd_pl;grd(lo:l)];
EE_pl = [EE_pl;EE(lo:l)];
EEE_pl = [real(EE_pl * conj(EE1/abs(EE1)) )-3*abs(delta),...
abs(EE_pl), real(EE_pl)-6*abs(delta),...
imag(EE_pl)-9*abs(delta)];  %-- plottting offset
EEE = [real(EE * conj(EE1/abs(EE1)) )-3*abs(delta),...
       abs(EE), real(EE)-6*abs(delta),...
       imag(EE)-9*abs(delta)];  %-- plotting offset
delr = real(delta);   deli = imag(delta);
axis('auto')
plot( grd_pl*2, EEE_pl, '-', [grd(1), 0.5]*2, abs(delta)*[1;1], '--',...
 [grd(1), 0.5]*2, abs(delta)*[1 -1;1 -1]-3*abs(delta), ':',...
 [grd(1), 0.5]*2, delr*[1 -1;1 -1]-6*abs(delta),':',...
 [grd(1), 0.5]*2, deli*[1 -1;1 -1]-9*abs(delta),':',...
 grd(jext)*2, EEE(jext,:), 'o', grd(iext)*2, EEE(iext,:), '+')
V = axis;
axis([-1*(~sym) 1 V(3) V(4)])
xlabel(getString(message('signal:crmz:NormalizedFrequency')))
title(getString(message('signal:crmz:WeightedError')))

%========================================================

function z_out = zeropad(x,num_of_zeros)
%ZEROPAD
%     zeropad(x,L) will append L zeros to each column of x.
%
[M,N] = size(x);
if M == 1,
   % input is a row vector
   z_out = [ x zeros(1,num_of_zeros) ];
else
   % input is a matrix or column
   z_out = [ x; zeros(num_of_zeros,N) ];
end

%========================================================

function y = db( x, dBrange, dBmax )
%DB Convert an array to decibels
%  Y=DB( X, dbRANGE, dbMAX ) will compute 20 Log(X)
%    and then scale or clip the result so that
%    the minimum dB level is dbMAX-dbRANGE.
%    ex: db(X, 80, 0) gives the range 0 to -80 dB
%  Y = dB( X, dbRANGE ) defaults to dBmax = 0
%
%  NOTE: on some machines log(0) gives an error
%        and aborts this function (esp. microVAX)

[M,N] = size(x);
if dBrange <= 0, error(message('signal:crmz:InvalidRange')); end
if nargin == 2, dBmax = 0; end
y = abs(x);
ymax = max(y)/10.0^(dBmax/20);
if M == 1,     %-- input is a row
   y = y/ymax;
else
   y = x ./ ymax(ones(M,1),:);
end 
thresh = 10.0^((dBmax-dBrange)/20);
y = y.*(y>thresh) + thresh.*(y<=thresh);
y = 20.0*log10(y);

%========================================================

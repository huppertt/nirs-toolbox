function [h,a,delta,HH,EE] = ...
     adesc( L, Lf, Lb, Ls, Lc, sym, Lfft, indx_edges, iext, HH, EE, a,...
            M_str, HH_str, h_str, A, delta, PLOTS, TRACE ) %#ok<*INUSL>
%ADESC Second-stage ascent-descent algorithm for CFIRPM
%  Solve problem on sets of extremal points using a descent method
%  based on the work of Demjanov & Malozemov and Wolfe. 

%   Authors: L. Karam, J. McClellan
%   Revised: 22-Oct-96, D. Orofino
%
%   Copyright 1988-2012 The MathWorks, Inc.

global DES_CRMZ WT_CRMZ GRID_CRMZ TGRID_CRMZ IFGRD_CRMZ %#ok<*NUSED>

% This parameter can be set to the desired accuracy:
ACCURACY = 0.01;

J = 1i;
is_odd = rem(L,2);
%Lfft = length(TGRID_CRMZ); % Passed

% Set the variable "bands".
vec_edges = TGRID_CRMZ(indx_edges);
bands = zeros(size(vec_edges));
for k = 1:length(vec_edges),
  if vec_edges(k)~=1,
    bands(k) = find(vec_edges(k)==GRID_CRMZ);
  else 
    bands(k) = length(GRID_CRMZ);
  end
end

no_stp = 1; 
a = a(1:Lb);
a = a(:);
n2 = length(a);
n = 2*length(a);
mxl = 2*L+1;
maj_it = 0;
alpha = 1;
nu = 0; 

% Set initial accuracy parameters:
acc_scale = max(abs(WT_CRMZ.*A)); 
acc_min =  ACCURACY;
acc = acc_min;
if acc < acc_min,
  acc = acc_min; 
end
r_min = 0.5;
r_o = acc_scale * r_min;
if r_o < r_min,
  r_o = r_min;
end
r = r_o;
epsi_min = ACCURACY;
epsi_o = acc_scale * epsi_min;
epsil = epsi_o; 
epsilon = 1-0.005;

% Construct initial subset:
HH_o = HH;
e_max = max(abs(EE)); 
[iext] = adesc_findset(EE,iext,bands,delta);
sub_EE = EE(iext);
sub_grd = GRID_CRMZ(iext);
sub_WT = WT_CRMZ(iext);
sub_max = max(abs(sub_EE));
fmax = iext;
[jext] = adesc_findextr(EE,bands,epsilon); 

% Check optimality of initial approximation:
gext = adesc_grad(EE, jext, GRID_CRMZ, WT_CRMZ, M_str, Lf, is_odd);
[rext] = adesc_minpolytope(gext, TRACE);
no_stp = (norm(rext) > acc);

% Begin major iteration:
while no_stp,

  jext     = find( (max(abs(sub_EE))-abs(sub_EE)) <= epsil );
  G        = adesc_grad(sub_EE, jext, sub_grd, sub_WT, M_str, Lf, is_odd);
  [NrG]    = adesc_minpolytope(G, TRACE);
  norm_NrG = norm(NrG);

  if (norm_NrG <= r)  
    nu = nu + 1; 
    epsil = epsi_o./(2^nu);
    r = r_o./(2^nu);
    f_extr = find(abs(sub_EE)/sub_max >= epsilon);
    gext = adesc_grad(sub_EE,f_extr,sub_grd,sub_WT,M_str,Lf,is_odd);
    [rext] = adesc_minpolytope(gext, TRACE);
    if norm(rext) < acc,
      maj_it = maj_it + 1;
      nu = 0;
      epsil = epsi_o; 
      r = r_o;  

      if TRACE,
        fprintf([getString(message('signal:adesc:Iter')) ': %3.0f   '],maj_it)
        fprintf([getString(message('signal:adesc:PeakWgtError')) ':  %5.2e  '],e_max)
        fprintf([getString(message('signal:adesc:SubsetPeakWgtErr')) ': %5.2e\n'],sub_max)
      end

      [iext] = adesc_findset(EE,iext,bands,sub_max);
      sub_EE = EE(iext);
      sub_grd = GRID_CRMZ(iext);
      sub_WT = WT_CRMZ(iext);
      sub_max = max(abs(sub_EE));

      if (length(iext) == length(fmax)),
        no_stp = ~all(iext==fmax);
      end

      fmax = iext;
    end

  else 
    d = - NrG./norm_NrG;
    d2 = d(1:n2)+J*d(n2+1:n);
    [HH_d] = adesc_reconst(d2,h_str,HH_str,M_str,TGRID_CRMZ,Lfft,...
                           L,Lf,Lc,Ls,Lb,is_odd,vec_edges,indx_edges );
    [HH_n,EE,e_max,alpha] = adesc_linsearch(alpha,iext,HH_o,HH_d,...
                                            sub_max,A,WT_CRMZ,IFGRD_CRMZ );
    a = a + alpha*d2;
    HH_o = HH_n;
    sub_EE = EE(iext);
    sub_max = max(abs(sub_EE));
  end

end % while no_stp

% end major iteration

HH = HH_o;
a = a(:); 
h = eval(h_str);
epsilon = 0.9*sub_max/e_max;
iext = adesc_findextr(EE,bands,epsilon);
jext = iext';
delta = max(abs(EE));

if TRACE,
  fprintf(getString(message('signal:adesc:OptimalSolutionObtained')))
end
% End Ascent-Descent Stage

%========================================================

function [fmax] = adesc_findextr(error,indx_edges,epsilon)
%ADESC_FINDEXTR
% returns indices of the extremals in "fmax".
%
error = error(:);
indx_edges = indx_edges(:);
Ngrid = length(error);
abs_e = abs(error);
abs_e = [abs_e(2); abs_e; abs_e(Ngrid-1)];
fmax = find( (abs_e(2:Ngrid+1) >= abs_e(1:Ngrid)) & ...
             (abs_e(2:Ngrid+1) > abs_e(3:Ngrid+2)) ); 
fmax = fmax(:);
fmax = sort([fmax; indx_edges]);
abs_e([1;Ngrid+2]) = [];
fmax(~(fmax(1:length(fmax)-1)-fmax(2:length(fmax))))=[];
if (length(fmax) > 1)
  fmax(abs_e(fmax)/max(abs_e) < epsilon) = []; 
end;

%========================================================

function [fmax] = adesc_findset(error,iext,indx_edges,delta)
% 
% Construct new set fmax of frequency points.  
% fmax: All local maxima including the set "iext".
%       Points with small deviation removed.
%
error = error(:);
indx_edges = indx_edges(:);
Ngrid = length(error);
abs_e = abs(error);
abs_e = [abs_e(2); abs_e; abs_e(Ngrid-1)];
fmax = find( (abs_e(2:Ngrid+1) >= abs_e(1:Ngrid)) & ...
             (abs_e(2:Ngrid+1) >= abs_e(3:Ngrid+2)) ); 
fmax = fmax(:);
fmax = sort([fmax; indx_edges; iext]);
abs_e([1;Ngrid+2]) = [];
fmax( fmax(1:length(fmax)-1) == fmax(2:length(fmax)) ) = [];
fmax( abs_e(fmax) < 0.9*abs(delta) ) = [];

%========================================================

function G = adesc_grad(EE, iext, grd, WT, M_str, Lf, is_odd)
%ADESC_GRAD Calculate gradients at the points iext.

J    = 1i;  % for use in M_str
fext = grd(iext);
W    = 2*pi*fext*((0:(Lf-1))+(~is_odd)*0.5);
Mb   = -eval(M_str)';
MWT  = diag(2.*WT(iext).*EE(iext));
M    = Mb*MWT;
G    = [real(M);imag(M)];

%========================================================

function [HHt,EEt,emxt,t] = ...
     adesc_linsearch(t0,iext,HH_o,d,emx,DD,WT,ifgrid)
%ADESC_LINSEARCH Performs a simple line search for step size. 
%  Efficiency of search can be improved by implementing 
%  method suggested in book "Introduction To Minimax" 
%  by V. F. Demjanov and V. N. Malozemov, pp. 109-112. 

t    = t0;
HH_o = HH_o(:);
d    = d(:);
no_stp = 1;
r_flg  = 0;
t_flg  = 1;
h_flg  = 0;
i_flg  = 0;
grtr   = 1;
c      = 2;
acc    = 1 - 10^(-0.1/10);
while grtr,
  HHt = HH_o + t*d; 
  EEt = WT .* (DD - HHt(ifgrid));
  emxt = max(abs(EEt(iext)));
  if (emxt <= emx),
    grtr = 0; 
  else
    t = t/2;
  end
end
tmin = t; HHt_min = HHt; EEt_min = EEt;
emin = emxt; I = 2*t*ones(2,1);
while no_stp,
  if (~i_flg), t = c*t; end
  HHt = HH_o + t*d; 
  EEt = WT .* (DD - HHt(ifgrid));
  emxt = max(abs(EEt(iext)));
  R = (emxt <= emin);
  if (R)
    tmin = t; HHt_min = HHt; EEt_min = EEt; emin = emxt;
    t_flg = 0; h_flg = 0;
    if (i_flg)
      I = sort([t;I(2)]);
      t = (t+I(2))/2;
    end
    r_flg = 1;
  elseif ((~R) && (i_flg))
      I = sort([I(1);t]);
      t = (t+I(1))/2;
  elseif (t_flg)
    t1 = t;
    t_flg = 0;
  elseif (r_flg)
    if (r_flg) 
      t1 = t;
    end;
    I = sort([tmin,t1]);
    t = (tmin+t1)/2;
    i_flg = 1;
  else
    no_stp = 0;
  end
  if  ( i_flg && ((I(2) - t) < acc*I(2)) ) 
    no_stp = 0;
  end;
end
HHt = HHt_min; EEt = EEt_min; emxt = max(abs(EEt)); t = tmin;

%========================================================

function [HH] = adesc_reconst(at,h_str,HH_str,M_str,tgrid,Lfft,... 
                L, Lf, Lc, Ls, Lb,is_odd, v_edges, in_edges)
%  Reconstructs the frequency response "HH" from the basis 
%  coefficients "at". 
%

global DES_CRMZ WT_CRMZ GRID_CRMZ TGRID_CRMZ IFGRD_CRMZ

J  = 1i;  % For possible use when evaluating HH_str
a  = at(:);
h  = eval(h_str);
hc = crmz_rotate(zeropad(h,Lfft-L), -Ls);
eval(HH_str);
W  = 2*pi*v_edges*((0:(Lf-1))+(~is_odd)*0.5);
Mb = eval(M_str);
HH(in_edges) = Mb * a;

%========================================================

function  Ptmin = adesc_minpolytope(P, TRACE)
%ADESC_MINPOLYTOPE
% Finds the nearest point Ptmin (smallest norm)
% in the convex hull of the set of points P1,...,PM
% given by the columns of the N x M matrix P.
%
% NOTE: the input matrix P is real-valued
%
% Implementation of a method by P. Wolfe presented in 
% "Mathematical Programming", vol. 11, 1976, pp. 128-149. 
%
[N,M] = size(P); 
if (M==1)
  Ptmin = P;
  return
end

Z1 = 1e-10;     % could be Z1 = 1e-12
Z2 = 1e-10;
Z3 = 1e-10;

P_norm = sum(abs(P).^2);
[pmn,J] = min(P_norm);
S = J;
w = 1;
no_stp = 1;
while (no_stp)
  X = P(:,S)*w;
  [pmn,J] = min(X'*P);
  PS_norm = max(sum(abs(P(:,S)).^2));
  PJ_norm = P(:,J)'*P(:,J);
  no_stp = (X'*P(:,J) <= (X'*X - Z1*max([PJ_norm;PS_norm])));
 if (no_stp)
    I = find(~(S-J), 1);
    if (~isempty(I))
      no_stp = 0;
      if TRACE, warning(message('signal:adesc:ComputingErr')); end
    end
 end
 if (no_stp)
  S = [S;J]; 
  w = [w;0];
  flg = 1;
  while (flg)
    e = ones(size(S));
    A = ones(length(S)) + P(:,S)'*P(:,S);
    u = A \ e;
    v = u./(e'*u);
    I = find(v <= Z2, 1);
    if isempty(I)
      w = v;
      flg = 0;
    else
      I = find((w-v) > Z3);
      t = w(I)./(w(I)-v(I));
      theta = min([1;1-min(t)]);
      w = theta*w + (1-theta)*v;
      I = find(w <= Z2);
      if (~isempty(I))
        w(I) = zeros(length(I),1);
        w(I(1)) = [];
        S(I(1)) = [];
      end
    end
  end
 end
end
Ptmin = X;

%========================================================

function z_out = zeropad(x,L)
%ZEROPAD ZEROPAD(X,L) will append L zeros to each column of X.
%
[M,N] = size(x);
if M==1,
   z_out = [x zeros(1,L)];
else
   z_out = [x; zeros(L,N)];
end

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
if M > 1              % ------- rotate columns ----------------
   num_places = mod(num_places,M);  % make num_places in range [0,M-1]
   rotated = [ x(M-num_places+1:M,:); x(1:M-num_places,:) ];
elseif N > 1          % ------- rotate row vector -------------
   num_places = mod(num_places,N);  % make num_places in range [0,N-1]
   rotated = [ x(N-num_places+1:N) x(1:N-num_places) ];
end

%========================================================

%#ok<*NASGU>
%#ok<*ASGLU>
%#ok<*AGROW>
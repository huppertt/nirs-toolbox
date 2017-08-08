function [ r, p, k ] = residuez( b, a, t )
%RESIDUEZ Z-transform partial-fraction expansion.
%   [R,P,K] = RESIDUEZ(B,A) finds the residues, poles and direct terms
%   of the partial-fraction expansion of B(z)/A(z),
%
%      B(z)       r(1)               r(n)
%      ---- = ------------ +...  ------------ + k(1) + k(2)z^(-1) ...
%      A(z)   1-p(1)z^(-1)       1-p(n)z^(-1)
%
%   B and A are the numerator and denominator polynomial coefficients,
%   respectively, in ascending powers of z^(-1).  R and P are column
%   vectors containing the residues and poles, respectively.  K contains
%   the direct terms in a row vector.  The number of poles is
%      n = length(A)-1 = length(R) = length(P)
%   The direct term coefficient vector is empty if length(B) < length(A);
%   otherwise,
%      length(K) = length(B)-length(A)+1
%
%   If P(j) = ... = P(j+m-1) is a pole of multiplicity m, then the
%   expansion includes terms of the form
%           R(j)              R(j+1)                      R(j+m-1)
%      -------------- + ------------------   + ... + ------------------
%      1 - P(j)z^(-1)   (1 - P(j)z^(-1))^2           (1 - P(j)z^(-1))^m
%
%   [B,A] = RESIDUEZ(R,P,K) converts the partial-fraction expansion back
%   to B/A form.
%
%   Warning: Numerically, the partial fraction expansion of a ratio of
%   polynomials represents an ill-posed problem.  If the denominator
%   polynomial, A(s), is near a polynomial with multiple roots, then small
%   changes in the data, including roundoff errors, can make arbitrarily
%   large changes in the resulting poles and residues. Problem
%   formulations making use of state-space or zero-pole representations
%   are preferable.
%
%   % Example:
%   %   Compute the partial fraction expansion of the following transfer
%   %   function H(z) = (1 + 2z^-1) / (1 - z^-1 + 2z^-2).
%
%   num = [1 1];                % Numerator coefficients
%   den = [1 -1 2];             % Denominator coefficients
%   [r,p] = residuez(num,den)   % H(z) = r(1)/(1-p(1)z^-1) + ...
%                               %        r(2)/(1-p(2)z^-1)
%
%   See also RESIDUE, PRONY, POLY, ROOTS, SS2TF, TF2SS, TF2ZP AND ZP2SS.

%   Author(s): J. McClellan, 10-24-90
%   Modified: T. Bryan, 10-22-97
%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] A.V. Oppenheim and R.W. Schafer, Discrete-Time
%         Signal Processing, Prentice-Hall, 1989, pp. 166-170.

mpoles_tol = 0.001;    %--- 0.1 percent
avg_roots = 0;         %--- FALSE, turns off averaging
if( nargin == 2 )
    %================= convert B/A to partial fractions =========
    if( a(1) == 0 )
        error(message('signal:residuez:SignalErr'))
    end
    % Cast to enforce Precision Rules
    if any([signal.internal.sigcheckfloattype(b,'single','residuez','B')...
        signal.internal.sigcheckfloattype(a,'single','residuez','A')])
      b = single(b);
      a = single(a);
    end
    
    b = b(:).'/a(1);   %--- b() & a() are now rows
    a = a(:).'/a(1);
    LB = length(b); 
    LA = length(a);
    if( LA == 1 )
        r = [];
        p = [];
        k = b;   %--- no poles!
        % Cast to enforce Precision Rules
        if isa(b,'single')
          r = single(r);
          p = single(p);
        end
        return
    end
    Nres = LA-1;
    % Cast to enforce Precision Rules
    p = zeros(Nres,1,class(b)); %#ok<*ZEROLIKE>
    r = p;
    LK = max([0;LB-Nres]);
    k = [];
    % Cast to enforce Precision Rules
    if isa(b,'single')
      k = single(k);
    end
    
    if( LK > 0 )
        [k,b] = deconv(fliplr(b),fliplr(a));
        k = fliplr(k);
        b(1:LK) = [];
        b = fliplr(b);
    end
    
    p = roots(a);
    N = Nres;         %--- number of terms in r()
    xtra = 1;         %--- extra pts invoke least-squares approx
    imp = zeros(N+xtra,1);     %--- zero padding to make impulse
    imp(1) = 1;
    h = filter( b, a, imp );
    % polmax = max(abs(p));
    polmax = 1.0;
    h = h.*filter(1,[1 -1/polmax],imp);  % if polmax == 1, this is just h = h;
    [mults, idx] = mpoles( p, mpoles_tol );
    p = p(idx);
    S = zeros(N+xtra,N);
    for j = 1:Nres
        if( mults(j) > 1 )            %--- repeated pole
            S(:,j) = filter( 1, [1 -p(j)/polmax], S(:,j-1) );
        else
            S(:,j) = filter( 1, [1 -p(j)/polmax], imp );
        end
    end
    jdone = zeros(Nres,1);
    bflip = fliplr(b);
    for i=1:Nres
        ii = 0;
        if( i==Nres )
            ii = Nres;
        elseif( mults(i+1)==1 )
            ii = i;
        end
        if( ii > 0 )
            if( avg_roots )       %---  average repeated roots ??
                jkl = (i-mults(i)+1):i;
                p(jkl) = ones(mults(i),1)*mean( p(jkl) );
            end
            temp = fliplr(poly( p([1:(i-mults(i)), (i+1):Nres]) ));
            r(i) = polyval(bflip,1/p(i)) ./ polyval(temp,1/p(i));
            jdone(i) = 1;
        end
    end
    %------ now do the repeated poles and the direct terms ------
    jjj = find(jdone==1);
    jkl = find(jdone~=1);
    if( any(jdone) )
        h = h - S(:,jjj)*r(jjj);
    end
    Nmultpoles = Nres - length(jjj);
    if( Nmultpoles > 0 )
        t = S(:,jkl)\h;
        r(jkl) = t(1:Nmultpoles);
    end
else
    %============= partial fractions ---> B(z)/A(z) ============
    %----- NOW, B <--> R
    %-----      A <--> P
    %-----      T <--> K
    % Cast to enforce Precision Rules
    if any([signal.internal.sigcheckfloattype(b,'single','residuez','R')...
        signal.internal.sigcheckfloattype(a,'single','residuez','P')...
        signal.internal.sigcheckfloattype(t,'single','residuez','K')])
      b = single(b);
      a = single(a);
      t = single(t);
    end
    
    res = b(:);     %--- r is now a column
    pol = a(:);     %--- p is now a column
    LR = length(res);
    LP = length(pol);
    LK = length(t);
    if( LR ~= LP )
        error(message('signal:residuez:InvalidDimensions'))
    end
    if( LP == 0 )
        p = []; r = t(:).';     %--- no poles!
        % Cast to enforce Precision Rules
        if isa(b,'single')
          p = single(p);
        end
        
        return
    end
    N = LP+LK;            %--- number of terms in b() and a()
    
    [mults, idx] = mpoles( pol, mpoles_tol, 0 );  % Maintain relative pole ordering
    pol = pol(idx);     %--- re-arrange poles & residues
    res = res(idx);
    p = poly(pol);  % p = p(:);          %--- p is really A(z)
    
    if( LK > 0 )
        r = conv( p, t );   %--- A(z)K(z)
    else
      % Cast to enforce Precision Rules
        r = zeros(1,N,class(b));     %--- r is B(z), returned as ROW
    end
    
    for i=1:LP
        temp = poly( pol([1:(i-mults(i)), (i+1):LP]) );
        r = r + [res(i)*temp zeros(1,N-length(temp))];
    end
    %  r=r(:).';
    %  p=p(:).';   % polynomials are ROWs
end

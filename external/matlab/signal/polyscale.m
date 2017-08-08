function y = polyscale(p,scale)
%POLYSCALE Scale roots of polynomial.
%   POLYSCALE(A,SCALE), where A is a vector of polynomial coefficients,
%   scales the roots of polynomial A by SCALE in the Z-plane.  If SCALE
%   is real and 0 <= SCALE <= 1, the roots of A are radially scaled
%   toward the origin in the Z-plane.  Complex SCALE provides arbitrary
%   changes to the root locations.
%
%   By reducing the radius of the roots in an autoregressive poly-
%   nomial, the bandwidth of the spectral peaks in the frequency
%   response is expanded (flattened).  This operation is often
%   referred to as "bandwidth expansion".
%
%   Example: Bandwidth expansion of LPC speech spectrum
%      load mtlb;                    % speech signal
%      Ao = lpc(mtlb(1000:1100),12); % 12th order AR polynomial
%      Ax = polyscale(Ao,.85);       % expanded bandwidth
%      subplot(2,2,1); zplane(1,Ao); title('Original');
%      subplot(2,2,3); zplane(1,Ax); title('Flattened');
%      [ho,w]=freqz(1,Ao);  [hx,w]=freqz(1,Ax);
%      subplot(1,2,2); plot(w,abs(ho), w,abs(hx));
%      legend('Original','Flattened');
%
%   See also POLYSTAB.

% Copyright 1988-2013 The MathWorks, Inc.

% The Z-transform P(z) is defined to be
%                    N-1
%   P(z) = Z{p(n)} = sum( p(n) z^(-n) )
%                    n=0
%
% Reducing the radius of the roots in the z-plane yields Px(z),
%   Px(z) = sum( p(n) * (scale*z)^(-n) )
%         = sum( (p(n) * scale^(-n)) * z^(-n) )
%         = sum( px(n) * z^(-n))
%         = Z{px(n)}
%
%  where px(n), the "bandwidth expanded" polynomial, is
%        px(n) = p(n) .* scale.^(-n)

% Cast to enforce precision rules
scale = signal.internal.sigcasttofloat(scale,'double','polyscale',...
  'SCALE','allownumeric');
% Checks if A is a valid numeric data input
signal.internal.sigcheckfloattype(p,'','polyscale','A');

y = p .* (scale .^ (0:length(p)-1));

% [EOF] polyscale.m

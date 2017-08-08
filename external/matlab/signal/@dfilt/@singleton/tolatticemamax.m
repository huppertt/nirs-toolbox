function Hd2 = tolatticemamax(Hd);
%TOLATTICEMAMAX  Convert to lattice MA maximum-phase.

%   Copyright 1988-2010 The MathWorks, Inc.
  
if ~isfir(Hd)
  error(message('signal:dfilt:singleton:tolatticemamax:DFILTErr'));
end
[b,a] = tf(Hd);
b = b/a(1);
if b(1)~=1,
    b = b/b(1);
    warning(message('signal:dfilt:singleton:tolatticemamax:GainIntroduced', num2str( 20*log10( 1/b( 1 ) ) )));
end
k = tf2latc(b,'max');
Hd2 = dfilt.latticemamax(k);

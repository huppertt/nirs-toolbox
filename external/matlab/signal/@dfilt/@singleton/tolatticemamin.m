
function Hd2 = tolatticemamin(Hd);
%TOLATTICEMAMIN  Convert to lattice moving-average.

%   Copyright 1988-2010 The MathWorks, Inc.
  
if ~isfir(Hd)
  error(message('signal:dfilt:singleton:tolatticemamin:DFILTErr'));
end
[b,a] = tf(Hd);
b = b/a(1);
if b(1)~=1,
    b = b/b(1);
    warning(message('signal:dfilt:singleton:tolatticemamin:GainIntroduced', num2str( 20*log10( 1/b( 1 ) ) )));
end
k = tf2latc(b);
Hd2 = dfilt.latticemamin(k);

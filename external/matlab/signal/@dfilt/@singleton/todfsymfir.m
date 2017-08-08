function Hd2 = todfsymfir(Hd);
%TODFSYMFIR  Convert to direct-form symmetric FIR.
%   Hd2 = TODFSYMFIR(Hd) converts discrete-time filter Hd to
%   direct-form symmetric FIR filter Hd2.

%   Copyright 1988-2002 The MathWorks, Inc.
  
if ~isfir(Hd)
  error(message('signal:dfilt:singleton:todfsymfir:DFILTErr'));
end
[b,a] = tf(Hd);
Hd2 = dfilt.dfsymfir(b/a);

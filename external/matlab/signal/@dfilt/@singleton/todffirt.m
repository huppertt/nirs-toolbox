function Hd2 = todffirt(Hd);
%TODFFIRT  Convert to direct-form FIR transposed.
%   Hd2 = (Hd) converts discrete-time filter Hd to direct-form FIR transposed
%   filter Hd2.

%   Copyright 1988-2002 The MathWorks, Inc.
  
if ~isfir(Hd)
  error(message('signal:dfilt:singleton:todffirt:DFILTErr'));
end
[b,a] = tf(Hd);
Hd2 = dfilt.dffirt(b/a(1));

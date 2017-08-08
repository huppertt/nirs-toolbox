function Hd2 = tolatticeallpass(Hd);
%TOLATTICEALLPASS  Convert to lattice allpass.
%   Hd2 = TOLATTICEALLPASS(Hd) converts discrete-time filter Hd to lattice
%   allpass filter Hd2.

%   Copyright 1988-2002 The MathWorks, Inc.
  
if ~isallpass(Hd),
    error(message('signal:dfilt:singleton:tolatticeallpass:DFILTErr'));
end

[b,a] = tf(Hd);
k = tf2latc(b,a);
Hd2 = dfilt.latticeallpass(k);

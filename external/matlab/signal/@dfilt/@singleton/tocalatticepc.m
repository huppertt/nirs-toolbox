function Hd2 = tocalatticepc(Hd);
%TOCALATTICEPC  Convert to couple allpass lattice, power complementary output.
%   Hd2 = TOCALATTICEPC(Hd) converts discrete-time filter Hd to couple
%   allpass lattice, power complementary filter Hd2.

%   Copyright 1988-2002 The MathWorks, Inc. 

[b,a] = tf(Hd);
if signalpolyutils('isfir',b,a),
    % FIR case
    if length(find(b~=0))>1,
        error(message('signal:dfilt:singleton:tocalatticepc:DFILTErr'));
    elseif length(find(b~=0))==1,
        % Pad a with zeros
        a(2:length(b)) = 0;
    end
end
if length(b)~=length(a),
    error(message('signal:dfilt:singleton:tocalatticepc:InvalidDimensions'));
end
[k1,k2,beta] = tf2cl(b,a);
Hd2 = dfilt.calatticepc(k1,k2,beta);

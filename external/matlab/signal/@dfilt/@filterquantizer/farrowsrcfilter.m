function [y,z,Tnext] = farrowsrcfilter(this,C,x,L,M,z,Tnext)
%FARROWSRCFILTER Filter implementation for MFILT.FARROWSRC

%   Copyright 2007 The MathWorks, Inc.

x = quantizeinput(this,x);
Lx = length(x);

    function updatestates(x)
        % Update the states for each scalar input x.
        z = [x;z(1:end-1,:)];
    end
y = [];
for k = 1:Lx,
    while Tnext>0
        y = [y;farrowfdfilter(this,C,x(k,:),Tnext/L,z)];
        Tnext  = Tnext - M;
    end
    updatestates(x(k,:));
    Tnext = Tnext + L;
end
end

% [EOF]

function C = dctmtx( t, fmax )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    N = size(t,1);
    Fs = 1/(t(2)-t(1));

    f = Fs/2/N * (0:N-1)'; % dct frq

    n = sum(f<fmax);
    [x,y] = ndgrid(0:N-1, 0:n-1);

    C = sqrt(2 / n) * cos(pi * (2*x + 1) .* y / (2 * n));
    C(:,1) = C(:,1) / sqrt(2);

end


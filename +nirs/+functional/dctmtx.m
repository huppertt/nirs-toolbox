function C = dctmtx( t, fmax )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    N = size(t,1);
    Fs = 1/(t(2)-t(1));

    f = Fs/2/N * (0:N-1)'; % dct frq

    n = sum(f <= fmax);
    
    [j,i] = ndgrid(0:n-1, 0:N-1);

    % DCT-II
    C = cos( pi/N * (i + 1/2) .* j );
    
    % normalization
    C = C*sqrt(2/N)';
    
    % orthogonalization
    C(1,:) = C(1,:) / sqrt(2);

end


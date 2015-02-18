function C = dctmtx( t, fmax )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    N = size(t,1);
    Fs = (t(2)-t(1));

    dct_t = (0:N-1)' / Fs;
    f = Fs/2/N * (0:N-1)'; % dct frq
    f = f( f < fmax );

    C = zeros(N,length(f));
    for i = 1:length(f);
        C(:,i) = cos(2*pi*f(i)*dct_t);
    end

end


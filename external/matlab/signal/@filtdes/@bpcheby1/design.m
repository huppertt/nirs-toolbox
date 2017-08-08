function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fpass1, Fpass2, Apass] = getdesignspecs(h, d);

if nargout == 1,
    hfdesign = fdesign.bandpass('N,Fp1,Fp2,Ap', d.Order, Fpass1, Fpass2, Apass);
    Hd       = cheby1(hfdesign);
    
    varargout = {Hd};
else

    % Set up design params
    N = get(d,'order');

    if rem(N,2),
        error(message('signal:filtdes:bpcheby1:design:MustBeEven'));
    end

    F = [Fpass1 Fpass2];

    [z,p,k] = cheby1(N/2,Apass,F);
    
    varargout = {z,p,k};
end

% [EOF]

function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fpass1, Fpass2, Apass, Astop] = getdesignspecs(h, d);

if nargout == 1,
    hfdesign = fdesign.bandstop('N,Fp1,Fp2,Ap,Ast', ...
        d.Order, Fpass1, Fpass2, Apass, Astop);
    Hd       = ellip(hfdesign);
    
    varargout = {Hd};
else
    % Set up design params
    N = get(d,'order');

    if rem(N,2),
        error(message('signal:filtdes:bsellip:design:MustBeEven'));
    end

    F = [Fpass1 Fpass2];

    [z,p,k] = ellip(N/2,Apass,Astop,F,'stop');
    varargout = {z,p,k};
end

% [EOF]

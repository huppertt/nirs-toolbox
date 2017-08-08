function varargout = design(h,d) %#ok<INUSL>
%DESIGN  Method to design the filter given the specs.

%   Copyright 1988-2012 The MathWorks, Inc.

if nargout == 1,
    hfdesign = fdesign.bandpass('N,F3dB1,F3dB2', d.Order, d.Fc1, d.Fc2);
    Hd       = butter(hfdesign);
        
    varargout = {Hd};
else
    % Set up design params
    N = get(d,'order');

    if rem(N,2),
        error(message('signal:filtdes:bpbutter:design:MustBeEven'));
    end

    % Get frequency specs, they have been prenormalized
    Fc1 = get(d,'Fc1');
    Fc2 = get(d,'Fc2');
    Fc = [Fc1, Fc2];

    [z,p,k] = butter(N/2,Fc);
    
    varargout = {z,p,k};

end

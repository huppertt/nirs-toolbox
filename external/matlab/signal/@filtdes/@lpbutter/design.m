function varargout = design(h,d) %#ok<INUSL>
%DESIGN  Method to design the filter given the specs.

%   Copyright 1988-2012 The MathWorks, Inc.

if nargout == 1,
    hfdesign = fdesign.lowpass('N,F3dB', d.Order, d.Fc);
    Hd       = butter(hfdesign);
    
    varargout = {Hd};
else
    [z,p,k] = butter(d.Order, d.Fc);
    
    varargout = {z,p,k};
end

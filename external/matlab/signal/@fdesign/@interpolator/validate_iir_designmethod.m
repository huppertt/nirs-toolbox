function validate_iir_designmethod(this,designMethod)
%VALIDATE_IIR_DESIGNMETHOD

%   Copyright 2009 The MathWorks, Inc.

if this.InterpolationFactor > 2
    response = this.Response;
    if strcmpi(response,'nyquist') || strcmpi(response,'halfband')
        if strcmpi(response(1),'h')
            response(1) = 'h';
        end        
        error(message('signal:fdesign:interpolator:validate_iir_designmethod:InvalidInterpFactor', designMethod, response));
    end
end

function varargout = validate_coeffs(this, coeffs)
%VALIDATE_COEFFS   Validate the coeffs

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

msg = [];

if length(coeffs) > 4 || length(coeffs) == 3,
    msgObj = message('signal:dfilt:wdfallpass:validate_coeffs:CoeffsMustBeVector');
    msg = getString(msgObj);
elseif length(coeffs) == 4 && (coeffs(1) ~= 0 || coeffs(3) ~= 0),
    msgObj = message('signal:dfilt:wdfallpass:validate_coeffs:InvalidValues');
    msg = getString(msgObj);
end
if isempty(msg),
    c = wdfcoefficients(this,coeffs);
    if any(abs(c.Section1) > 1) || any(isnan(c.Section1)),
      msgObj = message('signal:dfilt:wdfallpass:validate_coeffs:NotSupported');
      msg = getString(msgObj);      
    end
end

if nargout
    varargout = {msg};
else
    if ~isempty(msg), error(msgObj); end
end

% [EOF]

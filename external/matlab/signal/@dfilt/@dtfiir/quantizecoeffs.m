function quantizecoeffs(h,eventData)
% Quantize coefficients

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

% Quantize the coefficients

if isempty(h.refnum) || isempty(h.refden)
    return;
end

% Add check for a0 ~= 1
if ~strcmpi(class(h.filterquantizer), 'dfilt.filterquantizer'),
    if any(h.refden(1,1)~=1),
        error(message('signal:dfilt:dtfiir:quantizecoeffs:invalidA0'));
    end
end

[h.privnum,h.privden] = quantizecoeffs(h.filterquantizer,h.refnum,h.refden);

% [EOF]

function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2009 The MathWorks, Inc.

if ~isempty(S),

    % Check data type, quantize if needed
    S = validatestatesobj(Hd.filterquantizer, S);
    nsections = Hd.nsections;
    if rem(size(S.Numerator,2), nsections)~=0,
        error(message('signal:dfilt:df1tsos:thissetstates:InvalidNumeratorStateDimensions', nsections));
    end
    if rem(size(S.Denominator,2), nsections)~=0,
        error(message('signal:dfilt:df1tsos:thissetstates:InvalidDenominatorStateDimensions', nsections));
    end
    
    % Reshape to one column per channel	
	ns = 2*nsections;
    if ns>0,
        ncols = prod(size(S.Numerator))/ns; %#ok<PSIZE> numel doesn't work as expected in fixed point
        S.Numerator = reshape(S.Numerator,ns,ncols);
        S.Denominator = reshape(S.Denominator,ns,ncols);
    end
end

Hd.HiddenStates = S;


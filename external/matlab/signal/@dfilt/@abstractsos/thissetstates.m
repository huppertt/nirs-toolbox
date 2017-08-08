function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2009 The MathWorks, Inc.

if ~isempty(S),
    % Check data type, quantize if needed
    S = validatestates(Hd.filterquantizer, S);
    nsections = Hd.nsections;
    if rem(size(S,2), nsections)~=0,
        error(message('signal:dfilt:abstractsos:thissetstates:InvalidDimensions', nsections));
    end
	% Reshape to one column per channel	
	ns = 2*nsections;
    w=warning('off');
	ncols = prod(size(S))/ns; %#ok<PSIZE> numel doesn't work as expected in fixed point
    warning(w);
	S = reshape(S,ns,ncols);
end
Hd.hiddenstates = S;

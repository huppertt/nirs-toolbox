function [NL, PrevIPorts, PrevOPorts, mainparams]=statespacefootconnect(q,NL,H,mainparams);
%STATESPACEFOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

coeff = coefficients(H);
Amat = coeff{1};
dim = size(Amat,1);

% sum and gain
for m = 1 : (2*dim-1)*dim
    set(NL.nodes(m),'qparam','double');
end

inputports = zeros(1,dim);
for m = 1:dim
    inputports(m) = (2*dim-1)*dim + m;
end


outputports = zeros(1,dim);
for m = 1 : dim
    for n = 1 : 2 : 2*dim - 1
        tempind = (m - 1)*(2*dim-1) + n;
        if n == 1
            NL.connect(tempind,1,tempind+1,1);
        else
            NL.connect(tempind,1,tempind-1,2);
            if n < 2*dim - 1
                NL.connect(tempind-1,1,tempind+1,1);
            end
        end
        NL.connect(inputports((n+1)/2),1,tempind,1);
    end
    outputports(dim-m+1) = tempind - 1;
end


% specify the interstage connection
% last layer, therefore no next input and output 
PrevIPorts = [];
PrevOPorts = [];

for m = 1:dim
    PrevIPorts = [PrevIPorts, filtgraph.nodeport(inputports(m),1)];
    PrevOPorts = [PrevOPorts, filtgraph.nodeport(outputports(m),1)];
end


function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=statespacebodyconnect(q,NL,H,mainparams)
%STATESPACEBODYCONNECT specifies the connection and quantization parameters in the
%conceptual body stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

coeff = coefficients(H);
Bmat = coeff{2};
dim = size(Bmat,1);

inputind = 5;

PrevIPorts = [];
PrevOPorts = [filtgraph.nodeport(inputind,1)];
NextIPorts = [];
NextOPorts = [];


for m = 1 : dim
    tempind = (m-1)*5 + 1;
    set(NL.nodes(tempind),'qparam','double');
    NL.connect(inputind,1,tempind,1);
    NL.connect(tempind,1,tempind+1,1);
    
    tempind = (m-1)*5 + 2;
    set(NL.nodes(tempind),'qparam','double');
    NL.connect(tempind,1,tempind+1,1);
    NextIPorts = [NextIPorts filtgraph.nodeport(tempind,2)];
    
    tempind = (m-1)*5 + 3;
    NL.connect(tempind,1,tempind+1,1);
    NextOPorts = [NextOPorts filtgraph.nodeport(tempind,1)];
    
    
    tempind = (m-1)*5 + 4;
    set(NL.nodes(tempind),'qparam','double');
    if m == 1
        NL.connect(tempind,1,tempind+1+5,2);
    else
        NL.connect(tempind,1,tempind+1,1);
    end
    
    tempind = (m-1)*5 + 5;
    if m ~= 1
        set(NL.nodes(tempind),'qparam','double');
        if m ~= dim
            NL.connect(tempind,1,tempind+5,2);
        end
    end
end
PrevOPorts = [PrevOPorts filtgraph.nodeport(tempind,1)];


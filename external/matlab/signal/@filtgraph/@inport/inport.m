function inp = inport(Nindex,Sindex,From)
%inport Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(0,3,nargin,'struct'));

inp = filtgraph.inport;

if nargin > 0
    inp.nodeIndex = Nindex;
end

if nargin > 1 
    inp.selfIndex = Sindex; 
end

if nargin > 2
    inp.setfrom(From);
end

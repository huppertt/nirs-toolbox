function IP = copy(ip)
% copy method to force a deep copy.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));

IP = feval(str2func(class(ip)));

IP.nodeIndex = ip.nodeIndex;
IP.selfIndex = ip.selfIndex;

if ~isempty(ip.from)
    IP.from = copy(ip.from);
end 

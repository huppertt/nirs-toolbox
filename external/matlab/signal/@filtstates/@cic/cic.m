function h = cic(integrator,comb,nsections)
%CIC   Cascaded Integrator-Comb (CIC) filter states.
%   H = FILTSTATES.CIC constructs a default CIC filter states object.
%
%   H = FILTSTATES.DFIIR(INTEGRATORSTATES,COMBSTATES) constructs an object
%   and sets its 'INTEGRATOR' and 'COMB' properties to INTEGRATORSTATES and
%   COMB STATES respectively.  
%
%   Example #1, construct the default object
%   h = filtstates.cic
%
%   Example #2, construct an object with Integrator and Comb states
%   as vectors of zeros.
%   h = filtstates.cic(zeros(4,1),zeros(4,1));
%
%   See also FILTSTATES/INT.

%   Author(s): P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 3
    nsections = 1;
end

h = filtstates.cic;

if nargin
    if iscell(integrator)
        integrator = cell2obj(integrator, nsections);
    end

    h.Integrator = integrator;
    if nargin > 1
        if iscell(comb)
            comb = cell2obj(comb, nsections);
        end
        
        h.Comb = comb;
    end
end

% -------------------------------------------------------------------------
function hobj = cell2obj(celldata, nsections)

nchannels = length(celldata)/nsections;
celldata  = reshape(celldata, nsections, nchannels);
for indx = 1:nsections
    for jndx = 1:nchannels
        hobj(indx,jndx) = filtstates.state(celldata{indx,jndx});
    end
end

% [EOF]

function DGDF = scalardggen(q,Hd,coeffnames,doMapCoeffsToPorts)
%SCALARDGGEN Directed Graph generator for Discrete FIR

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(4,4,nargin,'struct'));

coefs = coefficients(reffilter(Hd));
num=coefs{1};

% Represent the filter in terms of DG_Dfilt
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;
DGDF = gen_DG_scalar_stages(q,num,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_dffir_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_scalar_stages(q,num,H,info,hTar)

info.nstages = 1;

Stg = scalarheader(q,num,H,info);

if info.doMapCoeffsToPorts
    Stg(2) = demux(q,H,1,info.coeffnames{1});
end

DGDF = filtgraph.dg_dfilt(Stg,'scalar');


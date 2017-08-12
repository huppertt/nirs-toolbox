function s = dataset2struct(d,varargin)
%DATASET2STRUCT Convert dataset array to structure array.
%   S = DATASET2STRUCT(D) converts the dataset array D to a structure array S.
%   Each variable of D becomes a field in S.  If D is an M-by-N array, then S
%   is M-by-1 and has N fields.  If D contains observation names, then S
%   contains those in the additional field 'ObsNames'.
%
%   S = DATASET2STRUCT(D,'AsScalar',true) converts the dataset array D to a
%   scalar structure S.  Each variable of D becomes a field in S.  If D is an
%   M-by-N array, the S has N fields, each of which has M rows.  If D contains
%   observation names, then S contains those in the additional field
%   'ObsNames'.
%
%   S = DATASET2STRUCT(D,'AsScalar',false) is identical to S = DATASET2STRUCT(D).
%
%   See also STRUCT2DATASET, DATASET2CELL, DATASET.

%   Copyright 2012 The MathWorks, Inc.


pnames = {'AsScalar'};
dflts =  {    false };
[asScalar] = statslib.internal.parseArgs(pnames, dflts, varargin{:});

[~,nvars] = size(d);
if asScalar
    s = cell2struct(d.data,d.varnames,2);
    if ~isempty(d.obsnames)
        s.ObsNames = d.obsnames;
        s = orderfields(s,[nvars+1 1:nvars]);
    end
else
    c = dataset2cell(d);
    s = cell2struct(c(2:end,:),c(1,:),2);
end

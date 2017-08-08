function [out coeffnames variables] = mapcoeffstoports(this,varargin)
%MAPCOEFFSTOPORTS 

%   Copyright 2009 The MathWorks, Inc.

out = parse_mapcoeffstoports(this,varargin{:});

stages = length(this.privallpasscoeffs); 

% default coefficient names
coeffnames = cell(1,stages);
for k = 1:stages
    coeffnames{k} = sprintf('C%d',k);
end

idx = find(strcmpi(varargin,'CoeffNames'));
if ~isempty(idx),
    userdefinednames = varargin{idx+1}; 
    % if user-defined coefficient names are empty, return the default names.
    if ~isempty(userdefinednames)
        coeffnames = userdefinednames;
    end
end

if length(coeffnames)~=length(this.privallpasscoeffs)
    error(message('signal:dfilt:cascadewdfallpass:mapcoeffstoports:InvalidParameter', length( this.privallpasscoeffs )));
end


% calculate actual gains
coefs = wdfcoefficients(this);
sectnames = fieldnames(coefs);
var = [];
for m = 1:length(sectnames)
    hap = coefs.(sprintf('Section%d',m));
    for n = 1:length(hap)
        gamma = hap(n);
        var(n) = local_wdfcoefficients(gamma);
    end
    variables{m} = var;
    var = [];
end

% If the coefficient variable of the first stage is empty, then there is
% no coefficient to be exported. Therefore, remove the coefficient name.
if isempty(variables{1})
    coeffnames = {};
end

%--------------------------------------------------------------------------
function coeffvar = local_wdfcoefficients(gamma)

coeffvar = 0;
if gamma ~= 0
    if gamma <= 1 && gamma > 0.5
        coeffvar = 1-gamma;
    elseif gamma <= 0.5 && gamma > 0
        coeffvar = gamma;
    elseif gamma < 0 && gamma >= -0.5
        coeffvar = -gamma;
    elseif gamma < -0.5 && gamma >= -1
        coeffvar = 1+gamma;
    end
end
% [EOF]

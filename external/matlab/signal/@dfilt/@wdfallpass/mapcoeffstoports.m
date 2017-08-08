function [out coeffnames variables] = mapcoeffstoports(this,varargin)
%MAPCOEFFSTOPORTS 

%   Copyright 2009 The MathWorks, Inc.

out = parse_mapcoeffstoports(this,varargin{:});

coeffnames = {'C'};
idx = find(strcmpi(varargin,'CoeffNames'));
if ~isempty(idx),
    userdefinednames = varargin{idx+1}; 
    % if user-defined coefficient names are empty, return the default names.
    if ~isempty(userdefinednames)
        coeffnames = userdefinednames;
    end
end

if length(coeffnames)~=1,
    error(message('signal:dfilt:wdfallpass:mapcoeffstoports:InvalidValue'));
end


% calculate actual gains
coefs = wdfcoefficients(this);
hap = coefs.('Section1');
var = zeros(1,length(hap));
for n = 1:length(hap)
    gamma = hap(n);
    var(n) = local_wdfcoefficients(gamma);
end

% If the coefficient is empty, do not export it to the workspace.
if isempty(var)
    variables = {};
    coeffnames = {};
else
    variables{1} = var;
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

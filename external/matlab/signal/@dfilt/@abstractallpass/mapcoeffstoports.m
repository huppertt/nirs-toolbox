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
    error(message('signal:dfilt:abstractallpass:mapcoeffstoports:InvalidValue'));
end

coeffs = this.AllpassCoefficients.';

% If the coefficient is empty, do not export it to the workspace.
if isempty(coeffs)
    variables = {};
    coeffnames = {};
else
    variables{1} = coeffs;
end

% [EOF]

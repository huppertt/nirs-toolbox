function [out coeffnames variables] = mapcoeffstoports(this,varargin)
%MAPCOEFFSTOPORTS 

%   Copyright 2008 The MathWorks, Inc.

out = parse_mapcoeffstoports(this,varargin{:});

coeffnames = {'Num'};
idx = find(strcmpi(varargin,'CoeffNames'));
if ~isempty(idx),
    userdefinednames = varargin{idx+1}; 
    % if user-defined coefficient names are empty, return the default names.
    if ~isempty(userdefinednames)
        coeffnames = userdefinednames;
    end
end

if length(coeffnames)~=1,
    error(message('signal:dfilt:dtffir:mapcoeffstoports:InvalidValue'));
end

Num = this.privNum.';
variables{1} = Num;

% [EOF]

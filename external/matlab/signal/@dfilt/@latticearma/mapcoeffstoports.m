function [out coeffnames variables] = mapcoeffstoports(this,varargin)
%MAPCOEFFSTOPORTS 

%   Copyright 2009 The MathWorks, Inc.

out = parse_mapcoeffstoports(this,varargin{:});

coeffnames = coefficientvariables(this);
idx = find(strcmpi(varargin,'CoeffNames'));
if ~isempty(idx), 
    userdefinednames = varargin{idx+1}; 
    % if user-defined coefficient names are empty, return the default names.
    if ~isempty(userdefinednames)
        coeffnames = userdefinednames;
    end
end

if length(coeffnames)~=2
    error(message('signal:dfilt:latticearma:mapcoeffstoports:InvalidValue'));
end

Lattice = this.privlattice.';
Ladder = this.privladder.';

% coefficients
variables{1} = Lattice;     
variables{2} = Ladder;
    
   

% [EOF]

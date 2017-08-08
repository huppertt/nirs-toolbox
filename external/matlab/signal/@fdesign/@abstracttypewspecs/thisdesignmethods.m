function [d, isfull, type] = thisdesignmethods(this, varargin)
%THISDESIGNMETHODS   Return the valid design methods.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

spec = get(this, 'CurrentSpecs'); 
if nargin > 1 
    if any(strcmpi(varargin{end}, set(this, 'Specification'))) 
        spec = feval(getconstructor(this, varargin{end})); 
        varargin(end) = []; 
    end 
end 

[d, isfull, type] = designmethods(spec, varargin{:});

% [EOF]

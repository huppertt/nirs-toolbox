function options = validateFinDiffRelStep(nVar,options,defaultopt)
%

%validateFinDiffRelStep sets the value for FinDiffRelStep.
%
%   Private to optimization solvers.

%   Copyright 2011 The MathWorks, Inc.

value = optimget(options,'FinDiffRelStep',defaultopt,'fast');

if isempty(value) % default is [] for most of the solvers
    % Check the finite difference type and set the appropriate value
    % Value of stepsize suggested in Trust Region Methods, Conn-Gould-Toint, section 8.4.3
    if strcmpi(options.FinDiffType,'forward')
        value = sqrt(eps)*ones(nVar,1); % Default for forward difference
    else
        value = eps^(1/3)*ones(nVar,1); % Default for central-diffference
    end
else
    value = value(:); % Make it a vector
    % User provided value; validate
    if ~isa(value,'double')
        ME = MException('optimlib:validateFinDiffRelStep:NotDoubleFinDiffRelStep', ...
            getString(message('optimlib:validateFinDiffRelStep:NotDoubleFinDiffRelStep')));
        throwAsCaller(ME);
    elseif isscalar(value) && value > 0
        % Valid scalar case; make it a vector
        value = value*ones(nVar,1);        
    elseif ~all(value > 0) % Must be positive
        ME = MException('optimlib:validateFinDiffRelStep:NotPositiveFinDiffRelStep', ...
            getString(message('optimlib:validateFinDiffRelStep:NotPositiveFinDiffRelStep')));
        throwAsCaller(ME);        
    elseif ~(length(value) == nVar)
        ME = MException('optimlib:validateFinDiffRelStep:InvalidSizeFinDiffRelStep', ...
            getString(message('optimlib:validateFinDiffRelStep:InvalidSizeFinDiffRelStep')));
        throwAsCaller(ME);
    end
    
end

% Set the value in the structure
options.FinDiffRelStep = value;

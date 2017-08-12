function options = createSolverOptions(solverName, varargin)
%

%CREATESOLVEROPTIONS Create a solver options object
%
%   OPTIONS = CREATESOLVEROPTIONS(SOLVERNAME, 'PARAM', VALUE, ...) creates
%   the specified solver options object. Any specified parameters are set
%   in the object.
%
%   See also OPTIMOPTIONS

%   Copyright 2012-2015 The MathWorks, Inc.

persistent optionsFactory 

% Create the optionsFactory map 
if isempty(optionsFactory) 
    % keys
    optimKeySet = {'linprog', ...
        'quadprog', ...
        'fminimax', ...
        'fmincon', ...
        'fgoalattain', ...
        'fseminf', ...
        'lsqnonlin', ...
        'lsqcurvefit', ...
        'fminunc', ...
        'fsolve', ...
        'lsqlin', ...
        'intlinprog'};
    globalOptimKeySet = {...
        'particleswarm', ...
        'pso'};
    % values
    optimValueSet = {@(varargin)optim.options.Linprog(varargin{:}), ...
        @(varargin)optim.options.Quadprog(varargin{:}), ...
        @(varargin)optim.options.Fminimax(varargin{:}), ...
        @(varargin)optim.options.Fmincon(varargin{:}), ...
        @(varargin)optim.options.Fgoalattain(varargin{:}), ...
        @(varargin)optim.options.Fseminf(varargin{:}), ...
        @(varargin)optim.options.Lsqnonlin(varargin{:}), ...
        @(varargin)optim.options.Lsqcurvefit(varargin{:}), ...
        @(varargin)optim.options.Fminunc(varargin{:}), ...
        @(varargin)optim.options.Fsolve(varargin{:}), ...
        @(varargin)optim.options.Lsqlin(varargin{:}), ...
        @(varargin)optim.options.Intlinprog(varargin{:})};
    globalOptimValueSet = {...
        @(varargin)optim.options.Particleswarm(varargin{:}), ...
        @(varargin)optim.options.Particleswarm(varargin{:})};
    
    % We check to see if there is an installation of Global Optimization
    % Toolbox here. Furthermore, we assume that these toolbox files will
    % not be removed between calls to this function. 
    %
    % Note that we do not perform a license check here to see if a user can
    % create a set of Global Optimization toolbox options. To ensure the
    % license check is correct, we have to check for license every time
    % this function is called. This is very expensive if optimoptions is
    % called multiple times in a tight loop. 
    %
    % As such, we rely on the license manager to throw an error in the case
    % where the user has a Global Optimization Toolbox installation, but
    % no license is available.
    hasGlobalOptimInstalled = ~isempty(ver('globaloptim'));
    
    % Create map
    if hasGlobalOptimInstalled
        optionsFactory = containers.Map(...
            [optimKeySet, globalOptimKeySet], ...
            [optimValueSet, globalOptimValueSet]);
    else
        optionsFactory = containers.Map(optimKeySet, optimValueSet);
    end
end

% Get creation function from factory
try
    optionsCreationFcn = optionsFactory(lower(solverName));
catch ME %#ok
    % handle unsupported solvers
    switch lower(solverName)
        case {'fminsearch', 'fzero', 'fminbnd', 'lsqnonneg'}
            upperSolver = upper(solverName);
            error(message('optimlib:options:createSolverOptions:MatlabSolversUnsupported', ...
                upperSolver, upperSolver));
        case 'ktrlink'
            ktrlink; % Call ktrlink to throw error message
        otherwise
            error(message('optimlib:options:createSolverOptions:InvalidSolver'));
    end
end

% Create the options
switch lower(solverName)
    case 'linprog'
        % Throw a more detailed error if a user tries to set old
        % options LargeScale or Simplex. Loop through the parameter names
        % and not the values.
        for i = 1:2:length(varargin)
            if ischar(varargin{i})
                if strcmpi(varargin{i}, 'LargeScale')
                    error(message('optimlib:options:createSolverOptions:LargeScaleUnsupported'));
                elseif strcmpi(varargin{i}, 'Simplex')
                    error(message('optimlib:options:createSolverOptions:SimplexUnsupported'));
                end
            end
        end
    case 'pso'
        error(message('optimlib:options:createSolverOptions:InvalidSolverName'));
end
% create options
options = optionsCreationFcn(varargin{:});

function options = optimoptions(solver, varargin)
%OPTIMOPTIONS  Create/modify optimization options
%
%   OPTIONS = OPTIMOPTIONS(SOLVER) creates optimization options, OPTIONS,
%   with the option parameters set to the default values relevant to the
%   optimization solver named in SOLVER, for example 'fmincon'. For a list
%   of valid optimization solvers, enter
%     doc optimoptions
%   to open the documentation page for optimoptions.
%
%   OPTIONS = OPTIMOPTIONS(SOLVER,'PARAM1',VALUE1,...) creates default
%   optimization options for SOLVER with the named parameters altered with
%   the specified values. OPTIONS can also be created using OPTIMOPTIONS
%   and setting the parameters via dot notation. For example:
%
%   % Create and set parameters using OPTIMOPTIONS
%   options = optimoptions('fminunc', 'TolX', 0.01, 'MaxIter', 1000)
%
%   % Create same options using OPTIMOPTIONS and dot notation
%   options = optimoptions('fminunc'); 
%   options.TolX = 0.01;
%   options.MaxIter = 1000;
%
%   OPTIONS = OPTIMOPTIONS(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of
%   OLDOPTS with the named parameters altered with the specified values.
%
%   OPTIONS = OPTIMOPTIONS(SOLVER, OTHEROPTS) first creates default
%   optimization options for SOLVER. The values of any parameters which are
%   common to OTHEROPTS and OPTIONS are copied to OPTIONS.
%
%   Note to see the parameters for a specific function, check the
%   documentation page for that function. For instance, enter
%     doc fmincon
%   to open the reference page for fmincon.
%
%      Example - setting options for a optimization solver.
%       % Create optimization options with the default parameters for
%       % FMINCON
%       options = optimoptions('fmincon')
%
%       % Set TolFun to 1e-3 and set MaxTime to 1000
%       options = optimoptions(options, 'TolFun', 1e-3, 'TolX', 1e-4)
%
%       % Set the Display option to 'iter'
%       options.Display = 'iter'
%
%       % Run the FMINCON solver with these options:
%       [x, fval, exitflag] = fmincon(@(x)peaks(x(1), x(2)), ...
%           [1 2], [1 1], 1, [], [], [-5 -5], [5 5], [], options)
%
%     Example - creating from another set of options
%
%      % Create some options for lsqnonlin
%      lsqOpts = optimoptions('lsqnonlin', 'MaxFunEvals', 1000)
%
%      % Use the above options for a call to lsqnonlin
%      x = lsqnonlin(@(x) sin(3*x),[1 4], [], [], lsqOpts)
%
%      % Now we want to add a nonlinear constraint to the above
%      % least squares problem. We will need to call fmincon, but we
%      % would like to use the same optimization options.
%
%      % First, create options for fmincon from the lsqnonlin options.
%      fmcOpts = optimoptions('fmincon', lsqOpts)
%     
%      % Call fmincon
%      [x, fval] = fmincon(@(x) (sin(3*x(1)))^2 + (sin(3*x(2)))^2, [1 4], ...
%             [], [], [], [], [], [], @(x)deal((x(1) - 1)^2 - x(2)), fmcOpts)

%   Copyright 2012-2014 The MathWorks, Inc.

% Error if no input is given
if nargin == 0
    error(message('optimlib:optimoptions:NotEnoughInputs'));
end

% Determine which solver we're creating options for. Note that invalid
% values of solver are caught in createSolverOptions.
if isa(solver, 'optim.options.SolverOptions')
    solverName = solver.SolverName;
    % If the user specified an old solver options object as the first
    % input, we need to pass this to the solver options object constructor.
    varargin = [{solver}, varargin]; 
elseif isa(solver, 'function_handle')
    solverName = func2str(solver);
else
    solverName = solver;
end

% Create the solver options object
try
    options = optim.options.createSolverOptions(solverName, varargin{:});
catch ME
    switch ME.identifier
        case 'MATLAB:InputParser:ParamMustBeChar'
            [linkToSolverCshStart, linkToSolverCshEnd] = ...
                i_createLinkToSolverCsh(solverName);
            error(message('optimlib:optimoptions:ParameterMustBeChar', ...
                upper(solverName), linkToSolverCshStart, linkToSolverCshEnd));
        case 'MATLAB:InputParser:UnmatchedParameter'
            [linkToSolverCshStart, linkToSolverCshEnd] = ...
                i_createLinkToSolverCsh(solverName);
            % Get the parameter that hasn't been matched
            idxQuote = regexp(ME.message, '''');
            error(message('optimlib:optimoptions:UnmatchedParameter', ...
                ME.message(idxQuote(1):idxQuote(2)), ...
                upper(solverName), upper(solverName), linkToSolverCshStart, linkToSolverCshEnd));
        case 'MATLAB:InputParser:ParamMissingValue'
            % Get the parameter that hasn't been matched
            idxQuote = regexp(ME.message, '''');
            error(message('optimlib:optimoptions:ParameterMissingValue', ...
                ME.message(idxQuote(1):idxQuote(2))));            
        otherwise
            throw(ME);
    end
end

function [docFileTagStart, docFileTagEnd] = i_createLinkToSolverCsh(solverName)
%i_createLinkToSolverCsh Create hyperlink tags to a solver doc page
%
%   [DOCFILETAGSTART, DOCFILETAGEND] = I_CREATELINKTOSOLVERCSH(SOLVERNAME)
%   creates hyperlink tags to the function reference documentation page for
%   the named solver. For example, for fmincon
%
%   docFileTagStart = <a href = "matlab: helpview('<DOCROOT>/toolbox/optim/ug/fmincon.html')">
%   docFileTagEnd = </a>
%
%   If links are not enabled in the MATLAB session, then the start and end
%   tags are empty strings.

% Don't add links to error/warning if running no-desktop, or deployed code, etc.
enableLinks = feature('hotlinks') && ~isdeployed;

% Create the hyperlink tags
if enableLinks
    if isempty(strfind(getfield(functions(str2func(solverName)),'file'),'globaloptim'))
        docFile = sprintf('%s/toolbox/optim/ug/%s.html', docroot, solverName);
    else
        docFile = sprintf('%s/toolbox/gads/%s.html', docroot, solverName);
    end
    docFileTagStart = sprintf('<a href = "matlab: helpview(''%s'')">', docFile);
    docFileTagEnd = sprintf('</a>');    
else
    docFileTagStart = '';
    docFileTagEnd = '';
end


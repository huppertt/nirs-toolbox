function [options, optionFeedback] = prepareOptionsForSolver(options, solverName)
%

%PREPAREOPTIONSFORSOLVER Prepare options for solver
%
%   [OPTIONS, OPTIONSFEEDBACK] = PREPAREOPTIONSFORSOLVER(OPTIONS,
%   SOLVERNAME) performs tasks to ensure the options are set up for use by
%   the solver. The following tasks are performed:-
%
%   * If a user has passed a SolverOptions object, first ensure that it is
%     a set of options for fseminf. 
%     
%   * If required by the solver, prepare strings to give feedback to users
%     on options they have or have not set. These are used in the exit
%     messages.
%
%   * If a user has passed a SolverOptions object, now extract the options
%     structure from the object.

%   Copyright 2012 The MathWorks, Inc.

% If a user has passed a SolverOptions object, first ensure that it is a
% set of options for the solver 
if isa(options, 'optim.options.SolverOptions')
    options = convertForSolver(options, solverName);
end

% If required, prepare strings to give feedback to users on options they
% have or have not set. These are used in the exit messages.
if nargout > 1
    optionFeedback = createOptionFeedback(options);
end

% If a user has passed a SolverOptions object, now extract the options
% structure from the object
if isa(options, 'optim.options.SolverOptions')
    options = extractOptionsStructure(options);
end

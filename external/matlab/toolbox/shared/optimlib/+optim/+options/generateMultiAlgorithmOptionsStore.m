function OS = generateMultiAlgorithmOptionsStore(OS)
%

%GENERATEMULTIALGORITHMOPTIONSSTORE Complete an OptionsStore structure for
%                                   a MultiAlgorithm options object
%
%   OS = GENERATEMULTIALGORITHMOPTIONSSTORE(OSIN) generates the fields of an
%   OptionsStore structure that can be automatically generated. OSIN must
%   contain the following fields:-
%
%   AlgorithmNames   : Cell array of algorithm names for the solver
%   DefaultAlgorithm : String containing the name of the default algorithm
%   AlgorithmDefaults: Cell array of structures. AlgorithmDefaults{i}
%                      holds a structure containing the defaults for 
%                      AlgorithmNames{i}.
%
%   GENERATEMULTIALGORITHMOPTIONSSTORE creates the following fields from
%   the above information:-
%
%   AlgorithmIndex      : 1-by-NumAlgorihms logical vector. The current
%                         algorithm is indicated by the position of true in
%                         the vector. The remaining elements are false.
%   IsConstantDefault   : Structure. Each field of the structure indicates
%                         whether the option has a default that varies for 
%                         different algorithms (false) or not (true).
%   SetByUser           : Structure. Indicate whether the user has set the 
%                         option (true) or not (false).
%   NonConstantDefaults : Structure. The fieldnames of the structure match 
%                         the names in the NonConstantDefaultFields field. 
%                         Each field contains a 1-by-NumAlgs cell array 
%                         which contains the default for that option for 
%                         each algorithm.
%   Options             : Options structure.

%   Copyright 2012 The MathWorks, Inc.

% Create the algorithm index
OS.AlgorithmIndex = strcmp(OS.DefaultAlgorithm, OS.AlgorithmNames);

% Number of algorithms
numAlgs = length(OS.AlgorithmDefaults);
% Create a list of all options
fnames = {};
for i = 1:numAlgs
    fnames = [fnames; fieldnames(OS.AlgorithmDefaults{i})]; %#ok
end
fnames = unique(fnames);

% Create the SetByUser structure
for i = 1:length(fnames)
    % SetByUser structure
    OS.SetByUser.(fnames{i}) = false;
end

% Create the NonConstantDefaults structure
OS = optim.options.createNonConstantDefaults(OS);

% Create the options structure
% First, set the values for the current algorithm options to their
% defaults.
OS.Options = OS.AlgorithmDefaults{OS.AlgorithmIndex};
% Next, find the options that do not apply to the current algorithm
RemOptions = setdiff(fnames, fieldnames(OS.Options));
% Next find the remaining algorithm indices
RemAlgs = find(~OS.AlgorithmIndex);
% For each of the remaining options 
for i = 1:length(RemOptions)
    if OS.IsConstantDefault.(RemOptions{i})
        for j = RemAlgs
            if isfield(OS.AlgorithmDefaults{j}, RemOptions{i})
                OS.Options.(RemOptions{i}) = OS.AlgorithmDefaults{j}.(RemOptions{i});
                break
            end
        end
    else
        OS.Options.(RemOptions{i}) = OS.NonConstantDefaults.(RemOptions{i}){OS.AlgorithmIndex};
    end
end

% Finally set the Algorithm option
OS.IsConstantDefault.Algorithm = true;
OS.SetByUser.Algorithm = false;
OS.Options.Algorithm = OS.DefaultAlgorithm;



classdef (Abstract) MultiAlgorithm < optim.options.SolverOptions
    %
    
    %MultiAlgorithm Options for Optimization Toolbox solvers with multiple
    %               algorithms
    %
    %   MultiAlgorithm is an abstract class representing a set of options
    %   for an Optimization Toolbox solver, where the solver has multiple
    %   algorithms. You cannot create instances of this class directly. You
    %   must create an instance of a derived class such
    %   optim.options.Fmincon.
    %
    %   See also OPTIM.OPTIONS.SOLVEROPTIONS
    
    %   Copyright 2012-2015 The MathWorks, Inc.
    
    properties (Dependent)
        
        %ALGORITHM Choose the optimization algorithm
        %
        %   For more information, see the "Options" section documentation page of
        %   the solver you are using.
        Algorithm
    end
    
    % Constructor
    methods (Hidden)
        function obj = MultiAlgorithm(varargin)
            %MULTIALGORITHM Construct a new MultiAlgorithm object
            %
            %   OPTS = OPTIM.OPTIONS.MULTIALGORITHM creates a set of options with each
            %   option set to its default value.
            %
            %   OPTS = OPTIM.OPTIONS.MULTIALGORITHM(PARAM, VAL, ...) creates a set of
            %   options with the named parameters altered with the specified values.
            %
            %   OPTS = OPTIM.OPTIONS.MULTIALGORITHM(OLDOPTS, PARAM, VAL, ...) creates
            %   a copy of OLDOPTS with the named parameters altered with the specified
            %   value
            
            % Call the superclass constructor
            obj = obj@optim.options.SolverOptions(varargin{:});
            
        end
    end
    
    
    methods
        function obj = set.Algorithm(obj, value)
            
            % Set Algorithm option
             obj = setProperty(obj, 'Algorithm', value, ...
                obj.OptionsStore.AlgorithmNames);            

            % If we get here, the property set has been successful and we
            % can update the OptionsStore
            
            % Set the algorithm index.
            obj.OptionsStore.AlgorithmIndex = strcmp(...
                obj.OptionsStore.AlgorithmNames, obj.OptionsStore.Options.Algorithm);            
            
            % Set the value of those options that have non-constant
            % defaults and haven't been set by the user
            for i = 1:length(obj.OptionsStore.NonConstantDefaultFields)
                if ~obj.OptionsStore.SetByUser.(obj.OptionsStore.NonConstantDefaultFields{i})
                    obj.OptionsStore.Options.(obj.OptionsStore.NonConstantDefaultFields{i}) = ...
                        obj.OptionsStore.NonConstantDefaults.(obj.OptionsStore.NonConstantDefaultFields{i}){obj.OptionsStore.AlgorithmIndex};
                end
            end
        end
        
        function value = get.Algorithm(obj)
            value = obj.OptionsStore.Options.Algorithm;
        end
        
    end
    
    % Display methods
    methods(Hidden, Access = protected)
        
        function header = addExtraHeader(obj, header)
            %ADDEXTRAHEADER Append text to the header
            %
            %   HEADER = ADDEXTRAHEADER(OBJ, HEADER) appends the following
            %   to the header in the following order:-
            %
            %     * The algorithm header
            %     * "Set by user" section of the display if no options are
            %     set by user.
            
            % Get algorithm header and append it to the main header
            algHeader = getAlgorithmHeader(obj);
            header = sprintf('%s%s', header, algHeader);
            
            % Call superclass addExtraHeader to append the extra header to
            % main header
            header = addExtraHeader@optim.options.SolverOptions(obj, header);
            
        end
        
        function footer = addExtraFooter(obj, footer)
            %ADDEXTRAFOOTER Add text to the footer
            %
            %   HEADER = ADDEXTRAFOOTER(OBJ, HEADER) prepends the following
            %   to the footer in the following order:-
            %
            %      * "Default" section of the display if all options are
            %      set by user.
            %      * The algorithm footer
            
            % Call superclass addExtraFooter to prepend the extra footer to
            % main footer
            footer = addExtraFooter@optim.options.SolverOptions(obj, footer);
            
            % Get algorithm footer and prepend it to the main footer
            algFooter = getAlgorithmFooter(obj);
            footer = sprintf('%s%s', footer, algFooter);
            
        end
        
        function algHeader = getAlgorithmHeader(obj)
            %GETALGORITHMHEADER Return the algorithm header for the display
            %
            %   HEADER = GETALGORITHMHEADER(OBJ) returns an algorithm specific string
            %   to be added to the header. An example of this string for fmincon is
            %   shown below:
            %
            %   "Options used by current Algorithm ('trust-region-reflective'):
            %    (Other available algorithms: 'active-set', 'interior-point', 'sqp')"
            
            
            % Now add the algorithm display header
            if isscalar(obj)
                
                %%% Create the first line of the algorithm display header
                currAlgStr = getString(message('optimlib:options:MultiAlgorithm:CurrentAlgorithmStr', ...
                    generateAlgorithmString, obj.Algorithm));
                algLine1 = sprintf('\n   %s', currAlgStr);
                
                %%% Create the second line of the algorithm display header
                
                % Create the "Other available algorithms" hyperlink
                if optim.options.SolverOptions.enableLinks
                    docFile = [docroot, '/optim/ug/', 'choosing-a-solver.html#bsbqd7i'];
                    docFileStartTag = sprintf('<a href = "matlab: helpview(''%s'')">', docFile);
                    docFileEndTag = '</a>';
                else
                    docFileStartTag = '';
                    docFileEndTag = '';
                end
                
                % Get a cell array of the other algorithm names
                otherAlgNames = obj.OptionsStore.AlgorithmNames(~obj.OptionsStore.AlgorithmIndex);
                
                % Put the other algorithm names in alphabetical order
                otherAlgNames = sort(otherAlgNames);
                
                % Create a comma-separated list of the other algorithm
                % names. We considered using cellfun to generate this list
                % but it is 3x slower than the following for loop. We also
                % looked at removing the if statement from the for loop and
                % adding the trailing brace afterwards. This was slightly
                % slower than the code below.
                otherAlgNamesStr = '';
                for i = 1:length(otherAlgNames)
                    if i < length(otherAlgNames)
                        otherAlgNamesStr = sprintf('%s''%s'', ', otherAlgNamesStr, otherAlgNames{i});
                    else
                        otherAlgNamesStr = sprintf('%s''%s''', otherAlgNamesStr, otherAlgNames{i});
                    end
                end
                
                % Get the "other available algorithms line of the header
                algLine2Text = getString(message('optimlib:options:MultiAlgorithm:OtherAvailableAlgorithmsStr', ...
                    docFileStartTag, docFileEndTag, otherAlgNamesStr));
                algLine2 = sprintf('   %s', algLine2Text);
                
                %%% Create the full algorithm header
                algHeader = sprintf('%s\n%s\n', algLine1, algLine2);
                
            else
                algHeader = '';
            end
            
        end
        
        function algFooter = getAlgorithmFooter(obj)
            %GETALGORITHMFOOTER Return the algorithm footer for the display
            %
            %   HEADER = GETALGORITHMHEADER(OBJ) returns an algorithm specific string
            %   to be added to the footer. An example of this string for fmincon is
            %   shown below:
            %
            %   "Show options not used by current Algorithm('trust-region-reflective')"
            
            % We're only customizing the display for scalar objects
            if isscalar(obj) && ~algorithmHasAllDisplayOptions(obj)
                
                % Create the algorithm footer
                
                % Below, we use evalc to capture the structure display for
                % the options that aren't used by the algorithm. evalc
                % always sets hotlinks to true (see g591312 for the
                % rationale for this). As such we need to capture the
                % current strong start and end tags and pass them to
                % showOptionsUnusedByAlgorithm.
                currentStrongStartTag = optim.options.SolverOptions.generateStrongStartTag; %#ok 
                currentStrongEndTag = optim.options.SolverOptions.generateStrongEndTag; %#ok
                
                % Generate the "Show options" link. Not keen on using evalc
                % here, but structs do not have a "toString" method. Note
                % that we hope to re-implement this in the future.
                showOptsStr = evalc('showOptionsUnusedByAlgorithm(obj, currentStrongStartTag, currentStrongEndTag)');
                if optim.options.SolverOptions.enableLinks
                    showOptsStr = regexprep(showOptsStr, '''', '''''');
                    showOptsStr = regexprep(showOptsStr, '\n', '\\n');
                    linkCmdStartTag = sprintf('<a href="matlab: fprintf(''%s'')">', ...
                        showOptsStr);
                    linkCmdEndTag = '</a>';
                    algFooterStr = getString(message('optimlib:options:MultiAlgorithm:ShowUnusedOptionsLinks', ...
                        linkCmdStartTag, linkCmdEndTag, ...
                        generateAlgorithmString, obj.Algorithm));
                    algFooter = sprintf('   %s\n', algFooterStr);
                else
                    algFooterStr = getString(message('optimlib:options:MultiAlgorithm:ShowUnusedOptionsNoLinks', 'Algorithm', obj.Algorithm));
                    algFooter = sprintf('   %s\n%s', algFooterStr, showOptsStr);
                end
                
                % Add a new line if all options are set by the user
                if needExtraFooter(obj)
                    algFooter = sprintf('\n%s', algFooter);
                end
                
            else
                algFooter = '';
            end
            
        end
        
        function allOptionsAtDefault = needExtraHeader(obj)
            %NEEDEXTRAHEADER Determine whether an extra header is needed
            %
            %   ALLOPTIONSATDEFAULT = NEEDEXTRAHEADER(OBJ) returns whether an extra
            %   header is needed.
            
            allOptionsAtDefault = (numCurrentAlgorithmDisplayOptionsSetByUser(obj) == 0);
        end
        
        function allOptionsSetByUser = needExtraFooter(obj)
            %NEEDEXTRAFOOTER Determine whether an extra footer is needed
            %
            %   ALLOPTIONSSETBYUSER = NEEDEXTRAFOOTER(OBJ) returns whether an extra
            %   footer is needed.
            
            allOptionsSetByUser = (numCurrentAlgorithmDisplayOptionsSetByUser(obj) == numCurrentAlgorithmDisplayOptions(obj));
        end
        
        function currentAlgorithmOptions = getAllDisplayOptions(obj)
            %GETDISPLAYOPTIONS Get the options to be displayed
            %
            %   CURRENTALGORITHMOPTIONS = GETDISPLAYOPTIONS(OBJ) returns a cell array of
            %   options to be displayed. For solver objects that inherit from
            %   MultiAlgorithm, options for the current algorithm are displayed.
            
            currentAlgorithmOptions = fieldnames(obj.OptionsStore.AlgorithmDefaults{obj.OptionsStore.AlgorithmIndex});
            currentAlgorithmOptions{end+1} = 'Algorithm';
        end
        
        
    end
    
    methods (Hidden)
        function OptionsStruct = extractOptionsStructure(obj)
            %EXTRACTOPTIONSSTRUCTURE Extract options structure from OptionsStore
            %
            %   OPTIONSSTRUCT = EXTRACTOPTIONSSTRUCTURE(OBJ) extracts a plain structure
            %   containing the options from obj.OptionsStore. The solver calls
            %   convertForSolver, which in turn calls this method to obtain a plain
            %   options structure.
            
            % Replace any special strings in the options object
            obj = replaceSpecialStrings(obj);
            
            % Extract main options structure
            OptionsStruct = obj.OptionsStore.Options;            
            
            % Optimoptions handles options with defaults that are algorithm
            % dependent. Some of these options may not be defined for the current
            % algorithm. If this is the case, we set the option to empty in the
            % options structure, which is the optimset standard.
            for i = 1:length(obj.OptionsStore.NonConstantDefaultFields)
                % Get the option name
                thisOption = obj.OptionsStore.NonConstantDefaultFields{i};
                % Is the default for thisOption defined for the current
                % algorithm? If it is not defined and the option has not
                % been set by the user, set the value to empty in the
                % structure.
                if ~isSetByUser(obj, thisOption) && ~isfield(obj.OptionsStore.AlgorithmDefaults{obj.OptionsStore.AlgorithmIndex}, thisOption)
                    OptionsStruct.(thisOption) = [];
                end
            end
            
            % Call method to map optimoptions for use in the solver and
            % optimtool
            OptionsStruct = mapOptionsForSolver(obj, OptionsStruct);            
            
        end
    end
    
    methods (Access = private)
        
        function numOptions = numCurrentAlgorithmDisplayOptions(obj)
            %NUMCURRENTALGORITHMDISPLAYOPTIONS Number of current algorithm
            %   options that are displayed
            %
            %   NUMOPTIONS = NUMCURRENTALGORITHMDISPLAYOPTIONS(OBJ) returns
            %   the number of options for the current algorithm that are
            %   displayed.
            
            allOptions = getDisplayOptions(obj);
            numOptions = length(allOptions);
        end
        
        function numSetByUser = numCurrentAlgorithmDisplayOptionsSetByUser(obj)
            %NUMCURRENTALGORITHMDISPLAYOPTIONSSETBYUSER Number of current
            %   algorithm display options that are set by the user
            %
            %   NUMSETBYUSER = NUMCURRENTALGORITHMDISPLAYOPTIONSSETBYUSER(OBJ) 
            %   returns the number of displayed options for the current 
            %   algorithm that are set by the user
            
            % Get names of all the properties for the current algorithm
            allOptions = getDisplayOptions(obj);
            numSetByUser = 0;
            for i = 1:length(allOptions)
                if obj.OptionsStore.SetByUser.(allOptions{i})
                    numSetByUser = numSetByUser + 1;
                end
            end
        end
        
        function hasAllOptions = algorithmHasAllDisplayOptions(obj)
%ALGORITHMHASALLDISPLAYOPTIONS Determine whether the algorithm has
%   all the options that the solver can display
%
%   HASALLOPTIONS = ALGORITHMHASALLDISPLAYOPTIONS(OBJ) returns true if the
%   set of options displayed for the current algorithm is identical to the
%   set of options that can be displayed by the solver.
            
            numAlgorithmOptions = numCurrentAlgorithmDisplayOptions(obj);
            numSolverOptions = length(properties(obj));
            hasAllOptions = (numAlgorithmOptions == numSolverOptions);
        end
    end
    
    methods (Access = private)
        
        function [setByUser, default] = getOptionsUnusedByAlgorithm(obj)
            %GETOPTIONSUNUSEDBYALGORITHM Return options unused by algorithm
            %
            %   [SETBYUSER, DEFAULT] = GETOPTIONSUNUSEDBYALGORITHM(OBJ) returns the
            %   options that are not used by the current algorithm. The cell array
            %   SETBYUSER returns the unused options that have been set by the user.
            %   The cell array, DEFAULT, returns the unused options that are at their
            %   default values.
            
            % Get the public (i.e. non-hidden) options for the solver
            allOptions = properties(obj);
            
            % Get the options for the current algorithm
            currAlgOptions = getDisplayOptions(obj);
            
            % Determine the options not used by the current algorithm
            allUnusedOptions = setdiff(allOptions, currAlgOptions);
            
            % Create the setByUser and default structures
            setByUser = struct([]);
            default = struct([]);
            for i = 1:length(allUnusedOptions)
                if obj.OptionsStore.SetByUser.(allUnusedOptions{i})
                    setByUser(1).(allUnusedOptions{i}) = obj.(allUnusedOptions{i});
                else
                    default(1).(allUnusedOptions{i}) = obj.(allUnusedOptions{i});
                end
            end
            
        end
    end
    
    methods
        function showOptionsUnusedByAlgorithm(obj, StrongStartTag, StrongEndTag)
            %SHOWOPTIONSUNUSEDBYALGORITHM Show options unused by algorithm
            %
            %   SHOWOPTIONSUNUSEDBYALGORITHM(OBJ, STRONGSTARTTAG,
            %   STRONGENDTAG) displays the options that are not used by the
            %   current algorithm.
            %
            %   See also OPTIMOPTIONS
            
            % Note that this method is called via evalc to capture the
            % structure display for the options that aren't used by the
            % algorithm. evalc always sets hotlinks to true (see g591312
            % for the rationale for this). As such this method requires the
            % tags from the calling workspace to be passed in, as these
            % tags reflect the hotlinks setting in the (actual) calling
            % workspace and not evalc's workspace.
                
            % Get the options unused by the current algorithm
            [setByUserOpts, defaultOpts] = getOptionsUnusedByAlgorithm(obj);
            
            % Display the unused options that have been set by the user
            if ~isempty(setByUserOpts)
                fprintf('   %s\n', ...
                    getString(message('optimlib:options:SolverOptions:SetByUserHeader', ...
                    StrongStartTag, StrongEndTag)));
                disp(setByUserOpts);
            end
            
            % Display the unused options that take their default value
            if ~isempty(defaultOpts)
                fprintf('   %s\n', getString(message('optimlib:options:SolverOptions:DefaultHeader', ...
                    StrongStartTag, StrongEndTag)));
                disp(defaultOpts);
            end
            
        end
        
    end
        
end


function algStr = generateAlgorithmString

StartTag = optim.options.SolverOptions.generateStrongStartTag;
EndTag = optim.options.SolverOptions.generateStrongEndTag;
algStr = sprintf('%sAlgorithm%s', StartTag, EndTag);
                
end

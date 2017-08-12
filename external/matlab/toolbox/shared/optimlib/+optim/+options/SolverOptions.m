classdef (Abstract) SolverOptions < matlab.mixin.CustomDisplay
    %
    
    %SolverOptions Base class for Optimization Toolbox options
    %   SolverOptions is an abstract class representing a set of options for an
    %   Optimization Toolbox solver. You cannot create instances of this class
    %   directly. You must create an instance of a derived class such
    %   optim.options.Fmincon.
    %
    %   For instances which are scalar, the object display is customized. We
    %   utilize the matlab.mixin.CustomDisplay class to help us create the
    %   display. To use this functionality, SolverOptions inherits from
    %   matlab.mixin.CustomDisplay.
    %
    %   See also MATLAB.MIXIN.CUSTOMDISPLAY
    
    %   Copyright 2012-2015 The MathWorks, Inc.
    
    properties (Hidden, Abstract, Access = protected)
        %OPTIONSSTORE Contains the option values and meta-data for the class
        %   OptionsStore is a structure containing a structure of option values for
        %   the solver. It also contains other meta-data that provides information
        %   about the state of the options, e.g. whether a given option has been
        %   set by the user.
        %
        %   This property must be defined and initialized by subclasses.
        OptionsStore
    end
    
    properties(Hidden, Access = protected)
        %VERSION Keep track of the version for the objects
        %   Version is an integer to indicate the version number of the object.
        Version
    end
  
    methods (Hidden, Access = protected)
        
        function obj = setProperty(obj, name, value, possValues)
            %SETPROPERTY Set a value for a solver option
            %
            %   OBJ = SETPROPERTY(OBJ, NAME, VALUE) sets the value of the named
            %   property to the specified value or throws an error.
            %
            %   OBJ = SETPROPERTY(OBJ, NAME, VALUE, POSSVALUES) allows callers to
            %   specify a list of possible values when the option is a cell string.
            
            if nargin < 4
                possValues = [];
            end
            if ischar(value)
                value = lower(deblank(value));
            end
            [value,validvalue, errmsg, errid] = ...
                optim.options.checkfield(name, value, possValues);
            if validvalue
                obj = setPropertyNoChecks(obj, name, value);
            else
                ME = MException(errid,'%s',errmsg);
                throwAsCaller(ME);
            end
        end
        
        function obj = setPropertyNoChecks(obj, name, value)
            %SETPROPERTYNOCHECKS Set a value for a solver option without
            %                    error checking
            %
            %   OBJ = SETPROPERTYNOCHECKS(OBJ, NAME, VALUE) sets the value
            %   of the named property to the specified value with no error
            %   checking. 
            % 
            %   NOTE: This method is designed for internal values of a
            %   public option only, e.g. setting the Display option to
            %   'testing'. For publicly supported values, use the
            %   setProperty method.
            
            obj.OptionsStore.Options.(name) = value;
            obj.OptionsStore.SetByUser.(name) = true;
            
        end
        
        function obj = upgradeFrom13a(obj, s)
            %UPGRADEFROM13A Perform common tasks to load 13a options
            %
            %   OBJ = UPGRADEFROM13A(OBJ, S) copies the OptionsStore and
            %   <Solver>Version fields from the structure into the object.
            %
            %   NOTE: Objects saved in R2013a will come in as structures.
            %   This is because we removed properties from
            %   optim.options.SolverOptions in 13b but did not need to
            %   write a loadobj at that time, as objects loaded correctly.
            %   Once we have a loadobj, the object loading throws an
            %   exception.

            % We need to reset OptionsStore
            obj.OptionsStore = s.OptionsStore;
                        
            % Set the version number back to allow loadobj to incrementally
            % upgrade the object correctly
            obj.Version = s.Version;
            
            % Note that s.SolverName is correctly set by the
            % constructor. The other three fields in s,
            % (StrongStartTag, StrongEndTag & EnableLinks) were removed
            % in R2013b.
        end
        
    end
    
    methods (Hidden)
        
        function obj = SolverOptions(varargin)
            %SOLVEROPTIONS Construct a new SolverOptions object
            %
            %   OPTS = OPTIM.OPTIONS.SOLVEROPTIONS creates a set of options with each
            %   option set to its default value.
            %
            %   OPTS = OPTIM.OPTIONS.SOLVEROPTIONS(PARAM, VAL, ...) creates a set of
            %   options with the named parameters altered with the specified values.
            %
            %   OPTS = OPTIM.OPTIONS.SOLVEROPTIONS(OLDOPTS, PARAM, VAL, ...) creates a
            %   copy of OLDOPTS with the named parameters altered with the specified
            %   value
            
            if nargin > 0
                
                % Deal with the case where the first input to the
                % constructor is a SolveOptions object.
                if isa(varargin{1},'optim.options.SolverOptions')
                    if strcmp(class(varargin{1}),class(obj))
                        obj = varargin{1};
                    else
                        % Get the properties from options object passed
                        % into the constructor.
                        thisProps = getOptionNames(varargin{1});
                        % Set the common properties. Note that we
                        % investigated first finding the properties that
                        % are common to both objects and just looping over
                        % those. We found that in most cases this was no
                        % quicker than just looping over the properties of
                        % the object passed in.
                        for i = 1:length(thisProps)
                            try %#ok
                                
                                % Store the value of option. Note that this
                                % is the default value because we haven't
                                % set this property yet.
                                defaultValue = obj.(thisProps{i});
                                
                                % Try to set one of the properties of the
                                % old object in the new one.
                                obj.(thisProps{i}) = varargin{1}.(thisProps{i});
                                
                                % If here, the property is common to both
                                % objects. We need to revert the SetByUser
                                % flag if the property has not been set and
                                % the default values of the property are
                                % equal.
                                %
                                % The following code is equivalent to this
                                % pseudo code:
                                % if isPropSetByUser
                                %      % If here, the property was set in
                                %      % the old object. In this case, as
                                %      % we have set the property in the
                                %      % new object there is nothing more
                                %      % to do here.
                                %  else
                                %      % If here, the property was not set
                                %      % in the old object. We now check
                                %      % whether the default value of the
                                %      % property in each object is
                                %      % identical.
                                %      if isDefaultEqual
                                %         % If the defaults are equal, then
                                %         % this property has not changed
                                %         % value. As such we do not want
                                %         % it to be marked as "SetByUser".
                                %      OptionsStore.SetByUser = false
                                %      else
                                %         % If the defaults are not equal
                                %         % we want to use the value from
                                %         % the old object. In this case,
                                %         % the property should be set and
                                %         % it should appear that it has
                                %         % been set. There is nothing more
                                %         % to do in this case.
                                %      end
                                % end
                                isPropSetByUser = ...
                                    varargin{1}.OptionsStore.SetByUser.(thisProps{i});
                                if ~isPropSetByUser && ...
                                        ( ischar(defaultValue) && strcmp(defaultValue, varargin{1}.(thisProps{i})) || ...
                                        isequal(defaultValue, varargin{1}.(thisProps{i})) )
                                    obj.OptionsStore.SetByUser.(thisProps{i}) = false;
                                end
                            end
                        end
                    end
                    firstInputObj = true;
                else
                    firstInputObj = false;
                end
                
                % Extract the options that the caller of the constructor
                % wants to set.
                if firstInputObj
                    pvPairs = varargin(2:end);
                else
                    pvPairs = varargin;
                end
                
                % Loop through each param-value pair and just try to set
                % the option. When the option has been fully specified with
                % the correct case, this is fast. The catch clause deals
                % with partial matches or errors.
                haveCreatedInputParser = false;
                for i = 1:2:length(pvPairs)
                    try
                        obj.(pvPairs{i}) = pvPairs{i+1};
                    catch ME %#ok
                        
                        % Create the input parser if we haven't already. We
                        % do it here to avoid creating it if possible, as
                        % it is slow to set up.
                        if ~haveCreatedInputParser
                            ip = inputParser;
                            % Structures are currently not supported as
                            % an input to optimoptions. Setting the
                            % StructExpand property of the input parser to
                            % false, forces the parser to treat the
                            % structure as a single input and not a set of
                            % param-value pairs.
                            ip.StructExpand =  false;
                            % Get list of option names
                            allOptionNames = getOptionNames(obj);
                            for j = 1:length(allOptionNames)
                                % Just specify an empty default as we already have the
                                % defaults in the options object.
                                ip.addParameter(allOptionNames{j}, []);
                            end
                            haveCreatedInputParser = true;
                        end
                        
                        % Get the p-v pair to parse.
                        thisPair = pvPairs(i:min(i+1, length(pvPairs)));
                        ip.parse(thisPair{:});
                        
                        % Determine the option that was specified in p-v pairs.
                        % These options will now be matched even if only partially
                        % specified (by 13a). Now set the specified value in the
                        % options object.
                        optionSet = setdiff(allOptionNames, ip.UsingDefaults);
                        obj.(optionSet{1}) = ip.Results.(optionSet{1});
                    end
                end
            end
        end
        
    end
    
    % Helper methods required for display
    methods (Hidden, Abstract, Access = protected)
        needExtraHeader(obj)
        needExtraFooter(obj)
        getAllDisplayOptions(obj)
    end
    
    % Display methods
    methods (Hidden, Access = protected)
        function dispOptions = getDisplayOptions(obj)
            %GETDISPLAYOPTIONS Get the options to be displayed
            %
            %   CURRENTALGORITHMOPTIONS = GETDISPLAYOPTIONS(OBJ) returns a cell array
            %   of options to be displayed. This function precludes the hidden options
            %   from the list.
            
            % First get the complete list of options
            dispOptions = getAllDisplayOptions(obj);
            
            % Now look for hidden or internal options and remove them from
            % the cell-array.
            
            % Query the metaclass to list the hidden options
            optsMetaInfo = metaclass(obj);
            nOpts = length(optsMetaInfo.PropertyList);
            
            % Go through each property, check for hidden properties
            % We also screen for dependent properties since we only care
            % about options
            % NOTE: we loop here since the list of properties is always >
            % in length than the options to be displayed
            for k = 1:nOpts
                if optsMetaInfo.PropertyList(k).Hidden && ...
                        optsMetaInfo.PropertyList(k).Dependent
                    dispOptions( ...
                        strcmp(optsMetaInfo.PropertyList(k).Name,dispOptions) ...
                        ) = [];
                end
            end
        end
        
        function header = addExtraHeader(obj, header)
            %ADDEXTRAHEADER Append text to the header
            %
            %   HEADER = ADDEXTRAHEADER(OBJ, HEADER) appends the "Set by user" section
            %   of the display to the header if no options are set by user.
            %
            %   There is currently no way of specifying a string to display when no
            %   options are set by the user. In this case, we have to just specify one
            %   property group and fold the Set by user group into the header.
            
            if needExtraHeader(obj)
                allDefaultStr = getString(message('optimlib:options:SolverOptions:AllDefaultStr'));
                header = sprintf('%s\n   %s\n     %s\n', ...
                    header, ...
                    getString(message('optimlib:options:SolverOptions:SetByUserHeader', ...
                    optim.options.SolverOptions.generateStrongStartTag, ...
                    optim.options.SolverOptions.generateStrongEndTag)), ...
                    allDefaultStr);
            end
            
        end
        
        function header = getHeader(obj)
            %GETHEADER Return the header for the display
            %
            %   HEADER = GETHEADER(OBJ) returns the header for the display. This method
            %   must be implemented as this class inherits from
            %   matlab.mixin.CustomDisplay.
            
            if isscalar(obj)
                if optim.options.SolverOptions.enableLinks
                    solverLink = sprintf('<a href="matlab: helpPopup %s" style="font-weight:bold">%s</a>', obj.SolverName, obj.SolverName);
                else
                    solverLink = obj.SolverName;
                end
                header = sprintf('  %s\n', ...
                    getString(message('optimlib:options:SolverOptions:HeaderStr', solverLink)));
                
                % Add extra header
                header = addExtraHeader(obj, header);
            else
                header = getHeader@matlab.mixin.CustomDisplay(obj);
            end
        end
        
        function footer = addExtraFooter(obj, footer)
            %ADDEXTRAFOOTER Add text to the footer
            %
            %   HEADER = ADDEXTRAFOOTER(OBJ, HEADER) prepends the "Default" section
            %   of the display to the footer if all options are set by user.
            %
            %   There is currently no way of specifying a string to display when no
            %   options are set by the user. In this case, we have to just specify one
            %   property group and fold the Set by user group into the header.
            
            if needExtraFooter(obj)
                allSetByUserStr = getString(message('optimlib:options:SolverOptions:AllSetByUserStr'));
                footer = sprintf('%s   %s\n     %s\n', ...
                    footer, ...
                    getString(message('optimlib:options:SolverOptions:DefaultHeader', ...
                    optim.options.SolverOptions.generateStrongStartTag, ....
                    optim.options.SolverOptions.generateStrongEndTag)), ...
                    allSetByUserStr);
            end
            
        end
        
        function footer = getFooter(obj)
            %GETFOOTER Return the footer for the display
            %
            %   FOOTER = GETFOOTER(OBJ) returns the footer for the display. This method
            %   must be implemented as this class inherits from
            %   matlab.mixin.CustomDisplay.
            
            if isscalar(obj)
                % Call superclass getFooter
                footer = getFooter@matlab.mixin.CustomDisplay(obj);
                
                % Add extra footer
                footer = addExtraFooter(obj, footer);
            else
                footer = getFooter@matlab.mixin.CustomDisplay(obj);
            end
        end
        
        function groups = getPropertyGroups(obj)
            %GETPROPERTYGROUPS Return the property groups for the display
            %
            %   GROUPS = GETPROPERTYGROUPS(OBJ) returns the property groups for the
            %   display. If all the options are at the default value, one group is
            %   returned containing all the options with the title "Default". If all
            %   the options have been set by the user, one group is returned with the
            %   title "Set by User". Otherwise two property groups are returned, one
            %   containing the options set by the user and another containing the
            %   remaining properties at their default values.
            %
            %   This method must be implemented as this class inherits from
            %   matlab.mixin.CustomDisplay.
            
            if isscalar(obj)
                
                % Get names of all properties to be displayed
                allOptions = getDisplayOptions(obj);
                
                % Sort the options to make them display in alphabetical
                % order.
                [~, idx] = sort(lower(allOptions));
                allOptions = allOptions(idx);
                
                % Create the group cell arrays
                if needExtraHeader(obj)
                    % Create the mixin property groups
                    defaultHeader = getString(message('optimlib:options:SolverOptions:DefaultHeader', ...
                        optim.options.SolverOptions.generateStrongStartTag, ...
                        optim.options.SolverOptions.generateStrongEndTag));
                    groups = matlab.mixin.util.PropertyGroup(allOptions, defaultHeader);
                elseif needExtraFooter(obj)
                    % Create the mixin property groups
                    setByUserHeader = getString(message('optimlib:options:SolverOptions:SetByUserHeader', ...
                        optim.options.SolverOptions.generateStrongStartTag, ...
                        optim.options.SolverOptions.generateStrongEndTag));
                    groups = matlab.mixin.util.PropertyGroup(allOptions, setByUserHeader);
                else
                    % Split the options into a group that have been set by
                    % user and a group taking the default value
                    idxSetByUser = false(1, length(allOptions));
                    for i = 1:length(allOptions)
                        idxSetByUser(i) = obj.OptionsStore.SetByUser.(allOptions{i});
                    end
                    setByUserGroup = allOptions(idxSetByUser);
                    defaultGroup = allOptions(~idxSetByUser);
                    
                    % Create the mixin property groups
                    setByUserHeader = getString(message('optimlib:options:SolverOptions:SetByUserHeader', ...
                        optim.options.SolverOptions.generateStrongStartTag, ...
                        optim.options.SolverOptions.generateStrongEndTag));
                    groups = matlab.mixin.util.PropertyGroup(setByUserGroup, setByUserHeader);
                    defaultHeader = getString(message('optimlib:options:SolverOptions:DefaultHeader', ...
                        optim.options.SolverOptions.generateStrongStartTag, ...
                        optim.options.SolverOptions.generateStrongEndTag));
                    groups(2) = matlab.mixin.util.PropertyGroup(defaultGroup, defaultHeader);
                end
            else
                groups = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            end
            
        end
        
    end
    
    % Required hidden helper methods
    methods (Abstract, Hidden)
        extractOptionsStructure(obj)
    end
    
    % Hidden helper methods
    methods (Hidden)
        function setByUserOptions = getOptionsSetByUser(obj)
            %GETOPTIONSSETBYUSER Return the options set by the user
            %
            %   SETBYUSEROPTIONS = GETOPTIONSSETBYUSER(OBJ) returns a cell string of
            %   option names that have been set by the user.
            
            allOptions = fieldnames(obj.OptionsStore.SetByUser);
            setByUser = struct2cell(obj.OptionsStore.SetByUser);
            setByUserOptions = allOptions([setByUser{:}]);
        end
        
        function obj = convertForSolver(obj, Solver)
            %CONVERTFORSOLVER Convert optimoptions for the named solver
            %
            %   OPTIONS = CONVERTFORSOLVER(OBJ, SOLVER) converts the supplied options
            %   to a set of options for the named solver.
            
            % It is possible for a user to pass in a vector of options to the
            % solver. Silently use the first element in this array.
            obj = obj(1);
            
            % Ensure Solver string is a proper name
            Solver(1) = upper(Solver(1));
            
            % Issue warning if user has passed an options object of the
            % wrong type.
            objStr = sprintf('optim.options.%s', Solver);
            if ~isa(obj, objStr)
                % Create the "Other available algorithms" hyperlink
                if optim.options.SolverOptions.enableLinks
                    docFile = [docroot, '/optim/ug/', 'optimoptions.html#btm5a5a-4'];
                    docFileStartTag = sprintf('<a href = "matlab: helpview(''%s'')">', docFile);
                    docFileEndTag = '</a>';
                else
                    docFileStartTag = '';
                    docFileEndTag = '';
                end
                warning(message('optimlib:options:SolverOptions:WrongSolverOptions', ...
                    upper(obj.SolverName), upper(Solver), ...
                    upper(Solver), upper(obj.SolverName), ...
                    docFileStartTag, upper(obj.SolverName), docFileEndTag));
                
                % Call the factory function to convert the solver object
                obj = optim.options.createSolverOptions(Solver, obj);
            end
            
        end
        
        function OptionsStruct = mapOptionsForSolver(~, OptionsStruct)
            %MAPOPTIONSFORSOLVER Map options for use by the solver
            %
            %   OptionsStruct = MAPOPTIONSFORSOLVER(obj, OptionsStruct)
            %   maps the specified structure so it can be used in the
            %   solver functions and in OPTIMTOOL.
            %
            %   This method does not alter the supplied structure.
            %   Subclasses can optionally supply mappings.
            
        end
        
        function [obj, OptimsetStruct] = mapOptimsetToOptions(obj, OptimsetStruct)
            %MAPOPTIMSETTOOPTIONS Map optimset structure to optimoptions
            %
            %   OBJ = MAPOPTIMSETTOOPTIONS(OBJ, OPTIMSETSTRUCT) maps specified optimset
            %   options, OptimsetStruct, to the equivalent options in the specified
            %   optimization object, obj.
            %
            %   [OBJ, OPTIONSSTRUCT] = MAPOPTIMSETTOOPTIONS(OBJ, OPTIMSETSTRUCT)
            %   additionally returns an options structure modified with any conversions
            %   that were performed on the options object.
            %
            %   This method does not alter the supplied options object or optimset
            %   structure. Subclasses can optionally supply mappings.
            
        end
        
        function OptionsStore = getOptionsStore(obj)
            %GETOPTIONSSTORE Return the OptionsStore
            %
            %   OPTIONSSTORE = GETOPTIONSSTORE(OBJ) returns the OptionsStore.
            
            OptionsStore = obj.OptionsStore;
        end
        
        function isSet = isSetByUser(obj, optionName)
            %ISSETBYUSER Return whether an option is set set by the user
            %
            %   ISSET = ISSETBYUSER(OBJ, OPTIONNAME) returns whether the specified
            %   option has been set by the user.
            
            isSet = obj.OptionsStore.SetByUser.(optionName);
        end
        
        function obj = replaceSpecialStrings(obj)
            %replaceSpecialStrings Replace special string values 
            %
            %   obj = replaceSpecialStrings(obj) replaces special string
            %   option values with their equivalent numerical value. We
            %   currently only use this method to convert FinDiffRelStep.
            %   However, in the future we would like to move the special
            %   string replacement code from the solver files to the
            %   options classes.
            %
            %   This method does not replace any special strings.
            %   Subclasses can optionally replace special strings.
            
        end
        
        function OptionNames = getOptionNames(obj)
            %GETOPTIONNAMES Get the options a user can set/get
            %
            %   OPTIONNAMES = GETOPTIONNAMES(OBJ) returns a list of the
            %   options that a user can set/get. This list is the union of
            %   the public properties and those that are hidden, dependent
            %   and have public set/get access.
            
            % All public properties are options that can be set/get by the
            % user.
            OptionNames = properties(obj);
            
            % Find any hidden options. These options will not be returned
            % from the call to properties. These options can still be
            % set/get by users.
            mc = metaclass(obj);
            numProps = length(mc.PropertyList);
            OptionsOnDeprecationPath = cell(numProps, 1);
            idxDep = false(1, numProps);
            for i = 1:length(mc.PropertyList)
                if mc.PropertyList(i).Dependent && mc.PropertyList(i).Hidden && ...
                        strcmp(mc.PropertyList(i).SetAccess, 'public') && ...
                        strcmp(mc.PropertyList(i).GetAccess, 'public')
                    OptionsOnDeprecationPath{i} = mc.PropertyList(i).Name;
                    idxDep(i) = true;
                end
            end
            OptionNames = [OptionNames; OptionsOnDeprecationPath(idxDep)];
            
        end
    end
    
    methods (Access=protected, Static)
        
        function StrongStartTag = generateStrongStartTag
        %GENERATESTRONGSTARTTAG Start tag for using bold face
        %
        %   STRONGSTARTTAG = GENERATESTRONGSTARTTAG returns a string
        %   holding bold face mark up ('<strong>') or no mark up ('')
        %   depending on whether we can use bold face mark up at the
        %   command prompt.
            
            if optim.options.SolverOptions.enableLinks
                StrongStartTag = '<strong>';
            else
                StrongStartTag = '';
            end
            
        end
        
        function StrongEndTag = generateStrongEndTag
        %GENERATESTRONGENDTAG End tag for using bold face
        %
        %   STRONGENDTAG = GENERATESTRONGENDTAG returns a string holding
        %   the end of the bold face mark up ('</strong>') or no mark up
        %   ('') depending on whether we can use bold face mark up at the
        %   command prompt.
            if optim.options.SolverOptions.enableLinks
                StrongEndTag = '</strong>';
            else
                StrongEndTag = '';
            end
        end
        
        function AreLinksEnabled = enableLinks
        %ENABLELINKS Indicate whether hyperlinks are enabled
        %
        %   ARELINKSENABLED = ENABLELINKS returns whether the customized
        %   display contains hyperlinks. EnableLinks holds a boolean flag
        %   indicating whether hyperlinks can be enabled or not.
        
        AreLinksEnabled = feature('hotlinks') && ~isdeployed;
        
        end
        
        
    end
   
end

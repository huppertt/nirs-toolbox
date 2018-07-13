classdef (Sealed) table
%TABLE Table.
%   Tables are used to collect heterogeneous data and metadata into a single
%   container.  Tables are suitable for storing column-oriented or tabular data
%   that are often stored as columns in a text file or in a spreadsheet.  Tables
%   can accommodate variables of different types, sizes, units, etc.  They are
%   often used to store experimental data, with rows representing different
%   observations and columns representing different measured variables.
%
%   Use the TABLE constructor to create a table from variables in the MATLAB
%   workspace.  Use the readtable function to create a table by reading data
%   from a text or spreadsheet file.
%
%   Tables can be subscripted using parentheses much like ordinary numeric
%   arrays, but in addition to numeric and logical indices, you can use a
%   table's variable and row names as indices.  You can access individual
%   variables in a table much like fields in a structure, using dot
%   subscripting.  You can access the contents of one or more variables using
%   brace subscripting.
%
%   Tables can contain different kinds of variables, including numeric, logical,
%   character, categorical, and cell.  However, a table is a different class
%   than the variables that it contains.  For example, even a table that
%   contains only variables that are double arrays cannot be operated on as if
%   it were itself a double array.  However, using dot subscripting, you can
%   operate on a variable in a table as if it were a workspace variable.  Using
%   brace subscripting, you can operate on one or more variables in a table as
%   if they were in a homogeneous array.
%
%   A table T has properties that store metadata such as its variable and row
%   names.  Access or assign to a property using P = T.Properties.PropName or
%   T.Properties.PropName = P, where PropName is one of the following:
%
%   TABLE metadata properties:
%       Description          - A string describing the table
%       DimensionNames       - A two-element cell array of strings containing names of
%                              the dimensions of the table
%       VariableNames        - A cell array containing names of the variables in the table
%       VariableDescriptions - A cell array of strings containing descriptions of the
%                              variables in the table
%       VariableUnits        - A cell array of strings containing units for the variables
%                              in table
%       RowNames             - A cell array of nonempty, distinct strings containing names
%                              of the rows in the table
%       UserData             - A variable containing any additional information associated
%                              with the table.  You can assign any value to this property.
%                       
%   TABLE methods and functions:
%     Construction and conversion:
%       table              - Create a table from workspace variables.
%       array2table        - Convert homogeneous array to table.
%       cell2table         - Convert cell array to table.
%       struct2table       - Convert structure array to table.
%       table2array        - Convert table to a homogeneous array.
%       table2cell         - Convert table to cell array.
%       table2struct       - Convert table to structure array.
%     Import and export:
%       readtable          - Create a table by reading from a file.
%       writetable         - Write a table to a file.
%       write              - Write a table to a file.
%     Size and shape:
%       istable            - True for tables.
%       size               - Size of a table.
%       width              - Number of variables in a table.
%       height             - Number of rows in a table.
%       ndims              - Number of dimensions of a table.
%       numel              - Number of elements in a table.
%       horzcat            - Horizontal concatenation for tables.
%       vertcat            - Vertical concatenation for tables.
%     Set membership:
%       intersect          - Find rows common to two tables.
%       ismember           - Find rows in one table that occur in another table.
%       setdiff            - Find rows that occur in one table but not in another.
%       setxor             - Find rows that occur in one or the other of two tables, but not both.
%       unique             - Find unique rows in a table.
%       union              - Find rows that occur in either of two tables.
%       join               - Merge two tables by matching up rows using key variables.
%       innerjoin          - Inner join between two tables.
%       outerjoin          - Outer join between two tables.
%     Data organization:
%       summary            - Print summary of a table.
%       sortrows           - Sort rows of a table.
%       stack              - Stack up data from multiple variables into a single variable.
%       unstack            - Unstack data from a single variable into multiple variables.
%       ismissing          - Find elements in a table that contain missing values.
%       standardizeMissing - Insert missing data indicators into a table.
%     Computations on tables:
%       varfun             - Apply a function to variables in a table.
%       rowfun             - Apply a function to rows of a table.
%
%   Examples:
%
%      % Create a table from individual workspace variables.
%      load patients
%      patients = table(LastName,Gender,Age,Height,Weight,Smoker,Systolic,Diastolic)
%
%      % Select the rows for patients who smoke, and a subset of the variables.
%      smokers = patients(patients.Smoker == true, {'LastName' 'Gender' 'Systolic' 'Diastolic'})
%
%      % Convert the two blood pressure variables into a single variable.
%      patients.BloodPressure = [patients.Systolic patients.Diastolic];
%      patients(:,{'Systolic' 'Diastolic'}) = []
%
%      % Pick out two specific patients by the LastName variable.
%      patients(ismember(patients.LastName,{'Smith' 'Jones'}), :)
%
%      % Convert the LastName variable into row names.
%      patients.Properties.RowNames = patients.LastName;
%      patients.LastName = []
%
%      % Use the row names to pick out two specific patients.
%      patients({'Smith' 'Jones'},:)
%
%      % Add metadata to the table.
%      patients.Properties.Description = 'Simulated patient data';
%      patients.Properties.VariableUnits =  {''  'Yrs'  'In'  'Lbs'  ''  'mm Hg'};
%      patients.Properties.VariableDescriptions{6} = 'Systolic/Diastolic';
%      summary(patients)
%
%      % Create a new variable in the table from existing variables.
%      patients.BMI = (patients.Weight * 0.453592) ./ (patients.Height * 0.0254).^2
%      patients.Properties.VariableUnits{'BMI'} =  'kg/m^2';
%      patients.Properties.VariableDescriptions{'BMI'} = 'Body Mass Index';
%
%      % Sort the table based on the new variable.
%      sortrows(patients,'BMI')
%
%      % Make a scatter plot of two of the table's variables.
%      plot(patients.Height,patients.Weight,'o')
%
%      % Create tables from text and spreadsheet files
%      patients2 = readtable('patients.dat','ReadRowNames',true)
%      patients3 = readtable('patients.xls','ReadRowNames',true)
%
%      % Create a table from a numeric matrix
%      load tetmesh.mat
%      t = array2table(X,'VariableNames',{'x' 'y' 'z'});
%      plot3(t.x,t.y,t.z,'.')
%
%   See also TABLE, CATEGORICAL

%   Copyright 2012-2014 The MathWorks, Inc.

     properties(Hidden = true,GetAccess='private', SetAccess='private')
        probe;
     end


    

    properties(Constant, GetAccess='private')
        version = 1.0;
        propsDflts = struct('Description'         , {''},...
                            'VariableDescriptions', {{}},...
                            'VariableUnits'       , {{}},...
                            'DimensionNames'      , {matlab.internal.table.dfltDimNames},...
                            'UserData'            , []);
        propertyNames = [ fieldnames(nirs.core.table.propsDflts); {'RowNames'; 'VariableNames'} ];
    end
    properties(GetAccess='private', SetAccess='private')
        ndims = 2;
        nrows = 0;
        rownames = {};
        nvars = 0;
        varnames = cell(1,0); % these can never be "truly" empty
        data = cell(1,0);
        
        % 'Properties' will also appear to contain 'VariableNames' and 'RowNames'.
        props =nirs.core.table.propsDflts;
    end
    properties(GetAccess='public', SetAccess='private', Dependent=true)
        %PROPERTIES Table metadata properties.
        %   T.Properties, where T is a table, is a scalar struct containing the
        %   following table metadata:
        %
        %       Description          - A string describing the table
        %       DimensionNames       - A two-element cell array of strings containing names of
        %                              the dimensions of the table
        %       VariableNames        - A cell array containing names of the variables in the table
        %       VariableDescriptions - A cell array of strings containing descriptions of the
        %                              variables in the table
        %       VariableUnits        - A cell array of strings containing units for the variables
        %                              in table
        %       RowNames             - A cell array of nonempty, distinct strings containing names
        %                              of the rows in the table
        %       UserData             - A variable containing any additional information associated
        %                              with the table.  You can assign any value to this property.
        %
        %   See also TABLE.
        Properties
    end
    methods
        function val = get.Properties(a), val = getProperties(a); end
    end

    methods
        function t = table(varargin)
%TABLE Create a table from workspace variables.
%   T = TABLE(VAR1, VAR2, ...) creates a table T from the workspace
%   variables VAR1, VAR2, ... .  All variables must have the same number
%   of rows.
%
%   T = TABLE(..., 'VariableNames', {'name1', ..., 'name_M'}) creates a
%   table containing variables that have the specified variable names.
%   The names must be valid MATLAB identifiers, and unique.
%
%   T = TABLE(..., 'RowNames', {'name1', ..., 'name_N'}) creates a table
%   that has the specified row names.  The names need not be valid MATLAB
%   identifiers, but must be unique.
%
%   Tables can contain variables that are built-in types, or objects that
%   are arrays and support standard MATLAB parenthesis indexing of the form
%   var(i,...), where i is a numeric or logical vector that corresponds to
%   rows of the variable.  In addition, the array must implement a SIZE method
%   with a DIM argument, and a VERTCAT method.
%
%
%   Examples:
%      % Create a table from individual workspace variables.
%      load patients
%      patients = table(LastName,Gender,Age,Height,Weight,Smoker,Systolic,Diastolic)
%      patients.Properties.Description = 'Simulated patient data';
%      patients.Properties.VariableUnits =  {''  ''  'Yrs'  'In'  'Lbs'  ''  'mm Hg'  'mm Hg'};
%
%   See also READTABLE, CELL2TABLE, ARRAY2TABLE, STRUCT2TABLE. 
        
        import matlab.internal.tableUtils.isstring
        
        if nargin == 0
            % nothing to do
            
        else
            % construct from separate variables
            vnames = repmat({''},1,nargin);
            argCnt = 0;

            % Set this as a first guess, it will shrink if there are name/value pairs.
            t_data = cell(1,nargin);

            varCnt = 0;
            numRows = 0;
            while argCnt < nargin
                argCnt = argCnt + 1;
                arg = varargin{argCnt};
                if isstring(arg) % matches any string, not just a parameter name
                    % Put that one back and start processing param name/value pairs
                    argCnt = argCnt - 1;
                    break
                elseif isa(arg,'function_handle')
                    error(message('MATLAB:table:FunAsVariable'));
                else % an array that will become a variable in t
                    varCnt = varCnt + 1;
                    t_data{varCnt} = arg;
                    vnames{varCnt} = inputname(argCnt); % might be empty
                end
                numRows_j = size(t_data{varCnt},1);
                if argCnt == 1
                    numRows = numRows_j;
                elseif ~isequal(numRows_j,numRows)
                    error(message('MATLAB:table:UnequalVarLengths'));
                end

            end % while argCnt < nargin, processing individual vars

            t.nrows = numRows;
            t.nvars = varCnt;
            t.data = t_data(1:varCnt);
            vnames = vnames(1:varCnt);

            if argCnt < nargin
                pnames = {'VariableNames'  'RowNames' };
                dflts =  {            []          []  };
                try
                [varnamesArg,rownamesArg,supplied] ...
                          = matlab.internal.table.parseArgs(pnames, dflts, varargin{argCnt+1:end});
                catch ME
                    % The inputs included a 1xM string that was interpreted as a param name,
                    % but something went wrong. If all of the preceedng inputs had one row,
                    % or if that string was first, these two errors suggest that the string
                    % was intended as data. Suggest alternative options, in that case.
                    if (t.nvars == 0 || t.nrows == 1)
                        if strcmp(ME.identifier,'MATLAB:table:parseArgs:BadParamName') || ...
                                strcmp(ME.identifier,'MATLAB:table:parseArgs:WrongNumberArgs')
                            cause = message('MATLAB:table:ConstructingFromStrings');
                            ME = addCause(ME,MException(cause.Identifier,'%s',getString(cause)));
                        end
                    end
                    throw(ME);
                end

                if supplied.VariableNames
                    vnames = varnamesArg;
                end
                if supplied.RowNames
                    if t.nvars == 0, t.nrows = length(rownamesArg); end
                    t = setRowNames(t,rownamesArg);
                end
            else
                supplied.VariableNames = false; supplied.RowNames = false;
            end

            if ~supplied.VariableNames
                % Fill in default names for data args where inputname couldn't
                empties = cellfun('isempty',vnames);
                if any(empties)
                    vnames(empties) = matlab.internal.table.dfltVarNames(find(empties)); %#ok<FNDSB>
                end
                % Make sure default names or names from inputname don't conflict
                vnames = matlab.lang.makeUniqueStrings(vnames,{},namelengthmax);
            end

            t = setVarNames(t,vnames); % error if invalid, duplicate, or empty
        end % if nargin > 0
        
        end % table constructor
    end % methods block
        
    
    methods(Hidden = true)
        % Variable Editor methods
        varargout  = variableEditorGridSize(t)
        [names,indices,classes,iscellstr,charArrayWidths] = variableEditorColumnNames(t)
        rowNames   = variableEditorRowNames(t)
        [code,msg] = variableEditorRowDeleteCode(t,workspaceVariableName,rowIntervals)
        [code,msg] = variableEditorColumnDeleteCode(t,workspaceVariableName,columnIntervals)
        t          = variableEditorPaste(t,rows,columns,data)
        t          = variableEditorInsert(t,orientation,row,col,data)
        [code,msg] = variableEditorSetDataCode(t,workspaceVariableName,row,col,rhs)
        [code,msg] = variableEditorUngroupCode(t,varName,col)
        [code,msg] = variableEditorGroupCode(t,varName,startCol,endCol)
        metaData   = variableEditorMetadata(t)
        [code,msg] = variableEditorMetadataCode(t,varName,index,propertyName,propertyString)
        [code,msg] = variableEditorRowNameCode(t,varName,index,rowName)
        [code,msg] = variableEditorSortCode(t,varName,tableVariableNames,direction)
        [code,msg] = variableEditorMoveColumn(t,varName,startCol,endCol)
        
        function d = getData(t)
            d = t.data;
        end
        
        
        
        
        
        % Allows tab completion after dot to suggest variables
        function p = properties(t)
            % This will be called for properties of an instance, but the
            % built-in will be still be called for the class name.  It will
            % return just Properties, which is correct.
            pp = [t.varnames(:); 'Properties'];
            if nargout == 0
                fprintf('%s\n',getString(message('MATLAB:ClassText:PROPERTIES_FUNCTION_LABEL','table')));
                fprintf('    %s\n',pp{:});
            else
                p = pp;
            end
        end
        function f = fieldnames(t), f = properties(t); end
        function f = fields(t),     f = properties(t); end
        
        
        
        function sz = numArgumentsFromSubscript(t,s,context)
            if isscalar(s) % one level of subscripting on a table
                sz = 1; % table returns one array for parens, braces, and dot
            elseif strcmp(s(end).type,'()')
                % This should never be called with parentheses as the last
                % subscript, but return 1 for that just in case
                sz = 1;
            else % multiple subscripting levels
                % Each level of subscripting should produce only one value, until
                % the last level. Get the value at the next to last level ...
                if strcmp(s(1).type,'{}')
                    x = subsrefBraces(t,s(1:end-1));
                elseif strcmp(s(1).type,'.')
                    x = subsrefDot(t,s(1:end-1));
                else % strcmp(s(1).type,'()'), e.g. t(...,...).Var
                    x = subsrefParens(t,s(1:end-1));
                end
                % ... and ask it how many values from the last level.
                if iscell(x)
                    sz = numel(x,s(end).subs{:});
                elseif isstruct(x)
                    sz = numel(x);
                elseif istable(x)
                    sz = 1; % table returns one array no matter what follows
                else % general case
                    try
                        sz = numArgumentsFromSubscript(x,s(end),context);
                    catch
                        if strcmp(s(end).type,'{}') % x{something}
                            sz = numel(x,s(end).subs{:});
                        else % strcmp(s(end).type,'.'), x.something, 
                            sz = numel(x);
                        end
                    end
                end
            end
        end
        
        % Methods we don't want to clutter things up with
        disp(t,bold,indent,fullChar)
        drawf= draw(t, vtype, vrange, thresh,figH)
        display(t)
        e = end(t,k,n)
        [varargout] = subsref(t,s)
        t = subsasgn(t,s,b)
        probe = getProbe(t)
        t = setProbe(t,probe)
        
        % Methods to override functions that we do not want to work
        function t = length(varargin),     throwUndefinedLengthError; end %#ok<STOUT>
        function t = transpose(varargin),  throwUndefinedError; end %#ok<STOUT>
        function t = ctranspose(varargin), throwUndefinedError; end %#ok<STOUT>
        function t = permute(varargin),    throwUndefinedError; end %#ok<STOUT>
        function t = reshape(varargin),    throwUndefinedError; end %#ok<STOUT>
        
        
        
    end % hidden methods block
    
    methods(Hidden = true, Static = true)
        function t = empty(varargin)
        %EMPTY Create an empty table.
        %   T = TABLE.EMPTY() creates a 0x0 table.
        %
        %   T = TABLE.EMPTY(NROWS,NVARS) or T = TABLE.EMPTY([NROWS NVARS]) creates an
        %   NROWSxNVARS table.  At least one of NROWS or NVARS must be zero.
        %
        %   See also TABLE, ISEMPTY.
            if nargin == 0
                t = table();
            else
                sizeOut = size(zeros(varargin{:}));
                if prod(sizeOut) ~= 0
                    error(message('MATLAB:table:empty:EmptyMustBeZero'));
                elseif length(sizeOut) > 2 && any(sizeOut(2:end) ~= 1)
                    error(message('MATLAB:table:empty:EmptyMustBeTwoDims'));
                else
                    % Create a 0x0 table, and then resize to the correct number
                    % of rows or variables.
                    t = table();
                    t.nrows = sizeOut(1);
                    t.nvars = sizeOut(2);
                    if t.nvars > 0
                        t.varnames = matlab.internal.table.dfltVarNames(1:t.nvars);
                        t.data = cell(1,t.nvars); % assume double
                    end
                end
            end
        end
        
        function t = loadobj(t)
            % Prevent fields in the props struct in future versions from getting through.
            toRemove = setdiff(fieldnames(t.props),fieldnames(table.propsDflts));
            if ~isempty(toRemove)
                t.props = rmfield(t.props,toRemove);
            end
        end
        
        % called by readtable
        t = readFromFile(filename,args)
        
        % called by cell2table, struct2table
        function t = fromScalarStruct(s)
            % Construct a table from a scalar struct
            vnames = fieldnames(s);
            p = length(vnames);
            if p > 0
                n = unique(structfun(@(f)size(f,1),s));
                if ~isscalar(n)
                    error(message('MATLAB:table:UnequalFieldLengths'));
                end
            else
                n = 0;
            end
            t = table;
            t.nrows = n;
            t.nvars = p;
            t.data = struct2cell(s)';
            t = setVarNames(t,vnames);
        end
    end % hidden static methods block
    
    methods(Access = 'protected')
        [varargout] = subsrefParens(t,s)
        [varargout] = subsrefBraces(t,s)
        [varargout] = subsrefDot(t,s)
        t = subsasgnParens(t,s,b,creating)
        t = subsasgnBraces(t,s,b)
        t = subsasgnDot(t,s,b)
        b = extractData(t,vars,like,a)
        t = replaceData(t,x,vars)
        c = extractRows(t,rowIndices,extractCells)
        t = replaceRows(t,rowIndices,c)
        
        writeTextFile(t,file,args)
        writeXLSFile(t,xlsfile,args)
        
        [group,glabels,glocs] = table2gidx(a,avars,reduce)
        
        function s = getProperties(t)
            %GET Get all table properties as a scalar struct.
            s = t.props; s.RowNames = t.rownames; s.VariableNames = t.varnames;
        end
        
        function t = setProperties(t,s)
            %SET Set some or all table properties from a scalar struct.
            if ~isstruct(s) || ~isscalar(s)
                error(message('MATLAB:table:InvalidPropertiesAssignment'));
            end
            fnames = fieldnames(s);
            for i = 1:length(fnames)
                fn = fnames{i};
                t = setProperty(t,fn,s.(fn));
            end
        end
        
        function b = cloneAsEmpty(a)
            %CLONEASEMPTY Create a new empty table from an existing one.
            if strcmp(class(a),'table') %#ok<STISA>
                b = table;
            else % b is a subclass of table
                b = a; % respect the subclass
                b.ndims = 2;
                b.nrows = 0;
                b.rownames = {};
                b.nvars = 0;
                b.varnames = cell(1,0); % these can never be "truly" empty
                b.data = cell(1,0);
                b.props = nirs.core.table.propsDflts;
            end
        end
    end % protected methods block
    
    methods(Static = true, Access = 'protected')
        t = readTextFile(file,args)
        t = readXLSFile(xlsfile,args)
        
        [ainds,binds] = table2midx(a,avars,b,bvars)
        [leftVars,rightVars,leftVarNames,rightVarNames,leftKeyVals,rightKeyVals,leftKeys,rightKeys,mergedKeys] ...
           = joinUtil(a,b,type,leftTableName,rightTableName,keys, ...
                      leftKeys,rightKeys,mergeKeys,leftVars,rightVars,ignoreDups,supplied)
        [c,il,ir] = joinInnerOuter(a,b,leftOuter,rightOuter,leftKeyvals,rightKeyvals, ...
                                   leftVars,rightVars,leftVarnames,rightVarnames)
    end % protected static methods block
end % classdef


%-----------------------------------------------------------------------------
function throwUndefinedLengthError
st = dbstack;
name = regexp(st(2).name,'\.','split');
m = message('MATLAB:table:UndefinedLengthFunction',name{2});
throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end % function

%-----------------------------------------------------------------------------
function throwUndefinedError
st = dbstack;
name = regexp(st(2).name,'\.','split');
m = message('MATLAB:table:UndefinedFunction',name{2});
throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end % function

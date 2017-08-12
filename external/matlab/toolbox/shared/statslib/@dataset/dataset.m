classdef dataset
%DATASET Dataset array.
%   Dataset arrays are used to collect heterogeneous data and metadata,
%   including variable and observation names, into a single container.  They
%   can be thought of as tables of values, with rows representing different
%   observations and columns representing different measured variables.
%   Dataset arrays are suitable for storing column-oriented or tabular data
%   that are often stored as columns in a text file or in a spreadsheet,
%   and can accommodate variables of different types, sizes, units, etc.
%
%   Use the DATASET constructor to create a dataset array from variables in
%   the MATLAB workspace.  You can also create a dataset array by reading data
%   from a text or spreadsheet file.  Dataset arrays can be subscripted using
%   parentheses much like ordinary numeric arrays, but in addition to numeric
%   and logical indices, you can use variable and observation names as
%   indices.  You can access each variable in a dataset array much like fields
%   in a structure, using dot subscripting.  Type "methods dataset" for a list
%   of operations available for dataset arrays.
%
%   Dataset arrays can contain different kinds of variables, including
%   numeric, logical, character, categorical, and cell.  However, a dataset
%   array is a different class than the variables that it contains.  For
%   example, even a dataset array that contains only variables that are double
%   arrays cannot be operated on as if it were itself a double array.
%   However, using dot subscripting, you can operate on variable in a dataset
%   array as if it were a workspace variable.
%
%   A dataset array D has properties that store metadata.  Access or assign to
%   a property using P = D.Properties.PropName or D.Properties.PropName = P,
%   where PropName is one of the following:
%
%       Description    - A string describing the data set
%       DimNames       - A two-element cell array of strings containing names of
%                        the dimensions of the data set
%       VarNames       - A cell array containing names of the variables in the data set
%       VarDescription - A cell array of strings containing descriptions of the variables
%                        in the data set
%       Units          - Units of variables in data set
%       ObsNames       - A cell array of nonempty, distinct strings containing names
%                        of the observations in the data set
%       UserData       - A variable containing additional information associated
%                        with the data set
%
%   Examples:
%      % Load a dataset array from a mat file and create some simple subsets
%      load hospital
%      h1 = hospital(1:10,:)
%      h2 = hospital(:,{'LastName' 'Age' 'Sex' 'Smoker'})
%
%      % Access and modify metadata
%      hospital.Properties.Description
%      hospital.Properties.VarNames{4} = 'Wgt'
%
%      % Create a new dataset variable from an existing one
%      hospital.AtRisk = hospital.Smoker | (hospital.Age > 40)
%
%      % Use individual variables to explore the data
%      boxplot(hospital.Age,hospital.Sex)
%      h3 = hospital(hospital.Age<30,{'LastName' 'Age' 'Sex' 'Smoker'})
%
%      % Sort the observations based on two variables
%      h4 = sortrows(hospital,{'Sex','Age'})
%
%   See also DATASET/DATASET, NOMINAL, ORDINAL

%   Copyright 2006-2013 The MathWorks, Inc.



    properties(Constant, GetAccess='private')
        propsFieldNames = ...
            {'Description'; 'VarDescription'; 'Units';   'DimNames'; 'UserData'};
        propsFieldDflts = ...
            {           '';               {};      {}; dfltDimNames;         []}
    end
    properties(GetAccess='private', SetAccess='private')
        ndims = 2;
        nobs = 0;
        obsnames = {};
        nvars = 0;
        varnames = cell(1,0); % these can never be "truly" empty
        data = cell(1,0);
        
        % 'Properties' will also appear to contain 'VarNames' and 'ObsNames'.
        props = cell2struct(dataset.propsFieldDflts, dataset.propsFieldNames, 1);
    end
    properties(GetAccess='public', SetAccess='private', Dependent=true)
        Properties;
    end
    methods
        function val = get.Properties(a), val = get(a); end
    end

    methods
        function a = dataset(varargin)
%DATASET Create a dataset array.
%   DS = DATASET(VAR1, VAR2, ...) creates a dataset array DS from the
%   workspace variables VAR1, VAR2, ... .  All variables must have the same
%   number of rows.
%
%   DS = DATASET(..., {VAR,'name'}, ...) creates a dataset variable named
%   'name' in DS.  Dataset variable names must be valid MATLAB identifiers,
%   and unique.
%
%   DS = DATASET(..., {VAR,'name1',...,'name_M'}, ...), where VAR is an
%   N-by-M-by-P-by-... array, creates M dataset variables in DS, each of size
%   N-by-P-by-..., with names 'name1', ..., 'name_M'.
%
%   DS = DATASET(..., 'VarNames', {'name1', ..., 'name_M'}) creates dataset
%   variables that have the specified variable names.  The names must be valid
%   MATLAB identifiers, and unique.  You may not provide both the 'VarNames'
%   parameter and names for individual variables.
%
%   DS = DATASET(..., 'ObsNames', {'name1', ..., 'name_N'}) creates a dataset
%   array that has the specified observation names.  The names need not be
%   valid MATLAB identifiers, but must be unique.
%
%   Dataset arrays can contain variables that are built-in types, or objects that
%   are arrays and support standard MATLAB parenthesis indexing of the form
%   var(i,...), where i is a numeric or logical vector that corresponds to
%   rows of the variable.  In addition, the array must implement a SIZE method
%   with a DIM argument, and a VERTCAT method.
%
%   You can also create a dataset array by reading from a text or spreadsheet
%   file, as described below.  This creates scalar-valued dataset variables,
%   i.e., one variable corresponding to each column in the file.  Variable
%   names are taken from the first row of the file.
%
%   DS = DATASET('File',FILENAME, ...) creates a dataset array by reading
%   column-oriented data in a tab-delimited text file.  The dataset variables
%   that are created are either double-valued, if the entire column is
%   numeric, or string-valued, i.e. a cell array of strings, if any element in
%   a column is not numeric.  Fields that are empty are converted to either
%   NaN (for a numeric variable) or the empty string (for a string-valued
%   variable).  Insignificant whitespace in the file is ignored.
%
%   Specify a delimiter character using the 'Delimiter' parameter name/value
%   pair.  The delimiter can be any of ' ', '\t', ',', ';', '|' or their
%   corresponding string names 'space', 'tab', 'comma', 'semi', or 'bar'.
%   Specify the number of lines to skip at the beginning of the file using the
%   'HeaderLines' parameter name/value pair.  Specify strings to be treated as
%   the empty string in a numeric column using the 'TreatAsEmpty' parameter
%   name/value pair.  This may be a character string, or a cell array of
%   strings.  'TreatAsEmpty' only applies to numeric columns in the file, and
%   numeric literals such as '-99' are not accepted.
%
%   DS = DATASET('File',FILENAME,'Format',FORMAT, ...)  creates a dataset
%   array using the TEXTSCAN function to read column-oriented data in a text
%   file.  Specifying a format can improve speed significantly for large
%   files.  FORMAT is a format string as accepted by the TEXTSCAN function.
%   You may also specify any of the parameter name/value pairs accepted by the
%   TEXTSCAN function, including the 'Delimiter' and 'HeaderLines' parameters.
%   The default delimiter when you specify the 'Format' parameter is ' '.
%
%   DS = DATASET('XLSFile',XLSFILENAME, ...) creates a dataset array from
%   column-oriented data in an Excel spreadsheet file.  You may also specify
%   the 'Sheet' and 'Range' parameter name/value pairs, with parameter values
%   as accepted by the XLSREAD function.  Variable names are taken from the
%   first row of the spreadsheet.  If the spreadsheet contains figures or
%   other non-tabular information, you should use the 'Range' parameter to
%   read only the tabular data.  By default, the 'XLSFile' option reads data
%   from the spreadsheet contiguously out to the right-most column that
%   contains data, including any empty columns that precede it.  If the
%   spreadsheet contains one or more empty columns between columns of data,
%   use the 'Range' parameter to specify a rectangular range of cells from
%   which to read variable names and data.
%
%   When reading from a text or spreadsheet file, the 'ReadVarNames' parameter
%   name/value pair determines whether or not the first row of the file is
%   treated as variable names.  Specify as a logical value (default true).
%
%   When reading from a text or spreadsheet file, the 'ReadObsNames' parameter
%   name/value pair determines whether or not the first column of the file is
%   treated as observation names.  Specify as a logical value (default false).
%   If the 'ReadVarNames' and 'ReadObsNames' parameter values are both true,
%   the name in the first column of the first row of the file is saved as the
%   first dimension name for the dataset.
%
%   DS = DATASET('XPTFile',XPTFILENAME, ...) creates a dataset array from a
%   SAS XPORT format file. Variable names from the XPORT format file are
%   preserved. Numeric data types in the XPORT format file are preserved
%   but all other data types are converted to cell arrays of strings. The
%   XPORT format allows for 28 missing data types. These are represented in
%   the file by an upper case letter, '.' or '_'. All missing data will be
%   converted to NaN values in DS. However, if you need the specific
%   missing types then you can recover this information using the XPTREAD
%   function. 
%
%   When reading from an XPORT format file, the 'ReadObsNames' parameter
%   name/value pair determines whether or not to try to use the first
%   variable in the file as observation names. Specify as a logical value
%   (default false). If the contents of the first variable are not valid
%   observation names then the variable will be read into a variable of the
%   dataset array and observation names will not be set.
%
%   Examples:
%      % Create a dataset array from workspace variables
%      load cereal
%      cereal = dataset(Calories,Protein,Fat,Sodium,Fiber,Carbo,Sugars, ...
%          'ObsNames',Name)
%      cereal.Properties.VarDescription = Variables(4:10,2);
%
%      % Create a dataset array from a single workspace variable
%      load cities
%      categories = cellstr(categories);
%      cities = dataset({ratings,categories{:}},'ObsNames',cellstr(names))
%
%      % Load data from a text or spreadsheet file
%      patients = dataset('File','hospital.dat','Delimiter',',','ReadObsNames',true)
%      patients2 = dataset('XLSFile','hospital.xls','ReadObsNames',true)
%
%   See also DATASET/SET, DATASET/GET, TDFREAD, TEXTSCAN,
%   XLSREAD, XPTREAD. 
        
        vnames = repmat({''},1,nargin);
        defaulted = true(1,nargin);
        argCnt = 0;
                
        if nargin == 0
            % nothing to do
            return
        elseif nargin == 1
            % convert a scalar struct
            arg = varargin{1};
            if isstruct(arg) && isscalar(arg)
                vnames = fieldnames(arg);
                p = length(vnames);
                if p > 0
                    n = unique(structfun(@(f)size(f,1),arg));
                    if ~isscalar(n)
                        error(message('stats:dataset:dataset:UnequalFieldLengths'));
                    end
                else
                    n = 0;
                end
                a.nobs = n;
                a.nvars = p;
                a.data = cell(1,p);
                for j = 1:p
                    a.data{j} = arg.(vnames{j});
                end
                a = setvarnames(a,vnames);
                return
            end
        end
        
        % Set these as a first guess, they may grow if some vars are
        % created by splitting up arrays, or shrink if there are
        % name/value pairs.
        a.nvars = nargin;
        a.data = cell(1,nargin);

        varCnt = 0;
        while argCnt < nargin
            argCnt = argCnt + 1;
            arg = varargin{argCnt};
            if isstring(arg)
                % Put that one back and start processing param name/value pairs
                argCnt = argCnt - 1;
                break
            elseif iscell(arg) && ~isscalar(arg) && isvector(arg) && (size(arg,1)==1)
                if (numel(arg)==2) && isstring(arg{2})
                    % {var,name}
                    if isa(arg{1},'dataset')
                        error(message('stats:dataset:dataset:DatasetVariable'));
                    end
                    varCnt = varCnt + 1;
                    a.data{varCnt} = arg{1};
                    vnames{varCnt} = arg{2};
                    defaulted(varCnt) = false;
                elseif (numel(arg)>2) && all(cellfun(@isstring,arg(2:end)))
                    % {var,name1,name2,...}
                    if isa(arg{1},'dataset')
                        error(message('stats:dataset:dataset:DatasetVariable'));
                    end
                    var = arg{1};                    
                    names = arg(2:end);
                    if length(names) ~= size(var,2)
                        error(message('stats:dataset:dataset:IncorrectNumberVarnames'));
                    end
                    szOut = size(var);
                    szOut(2) = []; if isscalar(szOut), szOut(2) = 1; end
                    ncols = size(var,2);
                    a.nvars = a.nvars + ncols-1;
                    vnames = [vnames repmat({''},1,ncols-1)]; %#ok<AGROW>
                    defaulted = [defaulted false(1,ncols-1)]; %#ok<AGROW>
                    a.data = [a.data cell(1,ncols-1)];
                    for j = 1:ncols
                        varCnt = varCnt + 1;
                        if ismatrix(var)
                            a.data{varCnt} = var(:,j);
                        else
                            a.data{varCnt} = reshape(var(:,j,:),szOut);
                        end
                        vnames{varCnt} = names{j};
                        defaulted(varCnt) = false;
                    end
                else
                    % false alarm -- cell-valued var without name
                    varCnt = varCnt + 1;
                    a.data{varCnt} = arg;
                    name = inputname(argCnt);
                    if isempty(name)
                        name = dfltvarnames(varCnt,true); % single name as char
                    end
                    vnames{varCnt} = name;
                end
            elseif isa(arg,'dataset')
                error(message('stats:dataset:dataset:DatasetVariable'));
            else
                % var without name
                varCnt = varCnt + 1;
                a.data{varCnt} = arg;
                name = inputname(argCnt);
                if isempty(name)
                    name = dfltvarnames(varCnt,true); % single name as char
                end
                vnames{varCnt} = name;
            end
            try
                nrows = size(a.data{varCnt},1);
            catch ME
                m = message('stats:dataset:dataset:SizeMethodFailed',varCnt);
                throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
            end
            if argCnt == 1
                a.nobs = nrows;
            elseif ~isequal(nrows,a.nobs)
                error(message('stats:dataset:dataset:UnequalVarLengths'));
            end

        end % while argCnt < nargin, processing individual vars
        
        a.nvars = varCnt;
        a.data = a.data(1:varCnt);
        vnames = vnames(1:varCnt);
        defaulted = defaulted(1:varCnt);
        % Make sure any names we filled in are not duplicates of non-empty names
        vnames = matlab.lang.makeUniqueStrings(vnames,defaulted,namelengthmax);
            
        if argCnt < nargin
            pnames = {'file' 'xlsfile' 'xptfile' 'varnames'  'obsnames' };
            dflts =  {    []        []        []         []          [] };
            [fileArg,xlsfileArg,xptfileArg,varnamesArg,obsnamesArg,supplied,otherArgs] ...
                      = dataset.parseArgs(pnames, dflts, varargin{argCnt+1:end});
            
            if supplied.file
                if argCnt > 0
                    error(message('stats:dataset:dataset:FileAndData'));
                end
                a = readFile(a,fileArg,otherArgs);
            elseif supplied.xlsfile
                if argCnt > 0
                    error(message('stats:dataset:dataset:XLSFileAndData'));
                end
                a = readXLSFile(a,xlsfileArg,otherArgs);
           elseif supplied.xptfile
                if argCnt > 0
                    error(message('stats:dataset:dataset:XPORTFileAndData'));
                end
                a = readXPTFile(a,xptfileArg,otherArgs);
            else
                if ~isempty(otherArgs)
                    error(message('stats:dataset:parseArgs:BadParamName',otherArgs{1}));
                end
            end
            
            if supplied.varnames
                if ~all(defaulted)
                    error(message('stats:dataset:dataset:VarNamesAndVarNamesParam'));
                else
                    vnames = varnamesArg;
                end
            end
            if supplied.obsnames
                if a.nvars == 0, a.nobs = length(obsnamesArg); end
                a = setobsnames(a,obsnamesArg);
            end
        end % if argCnt < nargin, processing name/value pairs
        
        % Varnames may be empty because we had no vars, or because we read
        % from a file.  In either case, no need to set them.
        if ~isempty(vnames)
            a = setvarnames(a,vnames); % names will be modified to make them valid
        end
        
        end % dataset constructor
    end % methods block
    
    methods(Hidden = true, Static = true, Access = 'public')
        function b = loadobj(b)
            % If loading an array from an old version without a VarDescription
            % property, fill in an empty value.
            if ~isfield(b.props,'VarDescription')
                b.props.VarDescription = {};
                % Put VarDescription second in the order
                fn = fieldnames(b.props);
                i = find(strcmp('VarDescription',fn));
                b.props = orderfields(b.props,[1 i 2:(i-1) (i+1):length(fn)]);
                
            % If loading an array that has a VarDescription property, but is
            % out of sync because the array was modified in an old version
            % that didn't know about VarDescription, clear out the property.
            else % isfield(b.props,'VarDescription')
                if numel(b.props.VarDescription) ~= b.nvars
                    b.props.VarDescription = {};
                end
            end
            
            % If loading an array that has a timeseries variable, need to
            % check if it is an "old-style" time series, and deal with it.
            % Regardless of its data length, an "old-style" object will have
            % been constructed as a 1x1 "new-style" array, and that probably
            % won't match the length of the dataset array.  Create a "new-style"
            % timeseries array sized to match the length of the dataset array,
            % leaving the original timeseries as the first element.  Pad out
            % with default elements, 
            isscalarTS = cellfun(@(c)isa(c,'timeseries') && isscalar(c),b.data,'UniformOutput',true);
            if any(isscalarTS)
                for i = find(isscalarTS)
                    ts = b.data{i};
                    if (ts.TimeInfo.Length == b.nobs) && (b.nobs ~= 1)
                        vn = b.varnames{i};
                        warning(message('stats:dataset:dataset:OldFormatTimeseries', vn, vn));
                        ts(b.nvars+1,1) = ts; ts = ts(1:end-1);
                        b.data{i} = ts;
                    end
                end
            end
            
            % Prevent fields in the props struct in future versions from getting through.
            toRemove = setdiff(fieldnames(b.props),dataset.propsFieldNames);
            if ~isempty(toRemove)
                b.props = rmfield(b.props,toRemove);
            end
        end
    end
    
    methods(Hidden = true)
        % Variable Editor methods
        varargout = variableEditorGridSize(a)
        [names,indices,classes,iscellstr,charArrayWidths] = variableEditorColumnNames(a)
        rowNames = variableEditorRowNames(a)
        [code,msg] = variableEditorRowDeleteCode(a,workspaceVariableName,rowIntervals)
        [code,msg] = variableEditorColumnDeleteCode(a,workspaceVariableName,columnIntervals)
        out = variableEditorPaste(a,rows,columns,data)
        out = variableEditorInsert(a,orientation,row,col,data)
        [code,msg] = variableEditorSetDataCode(a,workspaceVariableName,row,col,rhs)
        [code,msg] = variableEditorUngroupCode(this,varName,col)
        [code,msg] = variableEditorGroupCode(this,varName,startCol,endCol)
        [code,msg] = variableEditorMetadataCode(this,varName,index,propertyName,propertyString)
        [code,msg] = variableEditorRowNameCode(this,varName,index,obsName)
        [code,msg] = variableEditorSortCode(this,varName,datasetVariableNames,direction)
        [code,msg] = variableEditorMoveColumn(this,varName,startCol,endCol)
        out = variableEditorMetadata(this)
    
        function b = fieldnames(a)
            b = properties(a);
        end
        function b = properties(a)
            b = [a.varnames(:); properties('dataset')];
        end
        
        function sz = numArgumentsFromSubscript(~,~,~)
            % dataset returns one array for parens and dot, and is only prepared
            % to handle one cell at a time. For deeper subscripting, returning 1
            % is the same thing that numel has always done.
            sz = 1;
        end
        
        % Methods that we inherit, but do not want
        function a = transpose(varargin),  a = throwUndefinedError; end
        function a = ctranspose(varargin), a = throwUndefinedError; end
        function a = permute(varargin),    a = throwUndefinedError; end
        function a = reshape(varargin),    a = throwUndefinedError; end
        function a = fields(varargin),     a = throwUndefinedError; end
    end % hidden methods block
        
    methods(Static = true)
        function a = empty(varargin)
            if nargin == 0
                a = dataset;
            else
                sizeOut = size(zeros(varargin{:}));
                if prod(sizeOut) ~= 0
                    error(message('stats:dataset:empty:EmptyMustBeZero'));
                elseif length(sizeOut) > 2
                    error(message('stats:dataset:empty:EmptyMustBeTwoDims'));
                else
                    % Create a 0x0 dataset, and then resize to the correct number
                    % of observations or variables.
                    a = dataset();
                    a.nobs = sizeOut(1);
                    a.nvars = sizeOut(2);
                    if a.nvars > 0
                        a.varnames = dfltvarnames(1:a.nvars);
                        a.data = cell(1,a.nvars);
                    end
                end
            end
        end
    end % static methods block
    
    methods(Access = 'protected')
        [varargout] = subsrefParens(a,s)
        [varargout] = subsrefBraces(a,s)
        [varargout] = subsrefDot(a,s)
        a = subsasgnParens(a,s,b,creating)
        a = subsasgnBraces(a,s,b)
        a = subsasgnDot(a,s,b)
        b = extractData(a,vars)
        [ainds,binds] = dataset2idx(a,avars,b,bvars)
    end % private methods block
    
    methods(Access = 'private')
        a = readFile(a,file,args)
        a = readXLSFile(a,xlsfile,args)
        a = readXPTFile(a,xptfile,args)
        a = scalarRepmat(a,m,n)
    end % private methods block
    
    methods(Static = true, Access = 'private')
        [varargout] = parseArgs(pnames,dflts,varargin)
    end % static private methods block
            
end % classdef


%-----------------------------------------------------------------------------
function tf = isstring(s) % require a row of chars, or possibly ''
tf = ischar(s) && (isrow(s) || all(size(s)==0));
end % function


%-----------------------------------------------------------------------------
function throwUndefinedError
st = dbstack;
name = regexp(st(2).name,'\.','split');
m = message('stats:dataset:UndefinedFunction',name{2});
throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end % function


%-----------------------------------------------------------------------------
function names = dfltDimNames
names = { getString(message('stats:dataset:uistrings:DfltObsDimName')) ...
          getString(message('stats:dataset:uistrings:DfltVarDimName')) };
end % function

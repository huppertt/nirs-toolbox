classdef nnfcnInfo < handle

% Copyright 2010-2011 The MathWorks, Inc.
  
  properties (SetAccess = private)
    
    % Standard
    mfunction = '';
    name = '';
    title = '';
    description = '';
    type = '';
    typeName = '';
    version = 0.0;
    source = '';
    deprecated = false;
    parameters = [];
    defaultParam = nnetParam;
    
    % Hidden
    subfunctions = [];
    hasParameters = false;
    numParameters = 0;
    parameterDefaults = struct;
    parameterDefaultValues = {};
    parameterNames = {};
  end
  
  methods
    
    function x = nnfcnInfo(mfunction,name,type,version,subfunctions)
      
      % TODO - Accept parameters optionally too
      if (nargin < 5) || isempty(subfunctions)
        subfunctions = struct;
      end
      
      if ~ischar(mfunction) || (ndims(mfunction)>2) || (size(mfunction,1)~=1)
        error(message('nnet:FcnName:NotString'));
      end
      if ~ischar(name) || (ndims(mfunction)>2) || (size(mfunction,1)~=1)
        error(message('nnet:nnfcnInfo:Name'));
      end
      if ~ischar(type) || (ndims(mfunction)>2) || (size(mfunction,1)~=1)
        error(message('nnet:nnfcnInfo:Type'));
      end
      if ~isnumeric(version)||(numel(version)~=1)||~isfinite(version)||(version<=0)
        error(message('nnet:nnfcnInfo:Version'));
      end
      
      if strcmp(type,'nntype.type_fcn');
        % Avoid recursion
        type_name = 'Type Function';
      else
        typeInfo = feval(type,'info');
        type_name = typeInfo.name;
      end
      
      x.mfunction = mfunction;
      x.name = name;
      x.title = [name ' ' type_name];
      x.description = nnfcn.get_mhelp_title(mfunction);
      x.type = type;
      x.typeName = type_name;
      x.version = version;
      x.source = 'nnet'; % TODO - make an argument
      x.deprecated = false;
      
      x.subfunctions.mfunction = mfunction;
      x.subfunctions.self = str2func(x.mfunction);
      for i=fieldnames(subfunctions)'
        fn = i{1};
        x.subfunctions.(fn) = subfunctions.(fn);
      end
    end
    
    function disp(x)
      isLoose = strcmp(get(0,'FormatSpacing'),'loose');
      if (isLoose), fprintf('\n'), end
      disp([' Neural ' x.typeName ' Info'])
      if (isLoose), fprintf('\n'), end
      disp([nnlink.prop2link('mfunction') nnlink.fcn2strlink(x.mfunction)]);
      disp([nnlink.prop2link('type') nnlink.fcn2strlink(x.type)]);
      disp([nnlink.prop2link('name') nnstring.str2str(x.name)]);
      disp([nnlink.prop2link('typeName') nnstring.str2str(x.typeName)]);
      disp([nnlink.prop2link('title') nnstring.str2str(x.title)]);
      disp([nnlink.prop2link('description') nnstring.str2str(x.description)]);
      disp([nnlink.prop2link('version') num2str(x.version)]);
      disp([nnlink.prop2link('source') nnstring.str2str(x.source)]);
      disp([nnlink.prop2link('deprecated') nnstring.bool2str(x.deprecated)]);
      disp([nnlink.prop2link('subfunctions') nnlink.paramstruct2str(x.subfunctions,false)]);
      disp([nnlink.prop2link('parameters') nnlink.params2str(x.parameters)]);
      fprintf([nnlink.prop2link('defaultParam') nnlink.paramstruct2str(x.defaultParam) '\n']);
      if (isLoose), fprintf('\n'), end
    end
    
    function [s,err] = parameterStructure(x,values)
      if length(values) > x.numParameters
        error(message('nnet:Args:TooManyParam'));
      end
      values = [values x.parameterDefaultValues((length(values)+1):end)];
      s = cell2struct(values,x.parameterNames,2);
      err = nntest.param(x.parameters,s);
      if isempty(err)
        err = feval(x.mfunction,'check_param',s);
      end
      if ~isempty(err) && (nargout < 2),nnerr.throw('Param',err); end
    end
    
    function [s,err] = overrideStructure(x,s,values)
      if length(values) > x.numParameters
        error(message('nnet:Args:TooManyParam'));
      end
      for i=1:length(values)
        s.(x.parameterNames{i}) = values{i};
      end
      err = nntest.param(x.parameters,s);
      if isempty(err)
        err = feval(x.mfunction,'check_param',s);
      end
      if ~isempty(err) && (nargout < 2),nnerr.throw('Param',err); end
    end
    
    % NNET 6.0 Compatibility
    function pn = pnames(x)
      pn = x.parameterNames;
    end
    function pd = pdefaults(x)
      pd = x.defaultParam;
    end
    function pn = fpnames(x)
      pn = x.parameterNames;
    end
    function pd = fpdefaults(x)
      pd = x.defaultParam;
    end
    
  end
  
  % PROTECTED METHODS
  methods (Access = protected)
    
    function setupParameters(x,params)
      x.hasParameters = true;
      x.parameters = params;
      x.hasParameters = true;
      x.numParameters = length(x.parameters);
      x.parameterDefaults = struct;
      x.parameterDefaultValues = cell(1,x.numParameters);
      x.parameterNames = cell(x.numParameters,1);
      for i=1:x.numParameters
        p = x.parameters(i);
        x.parameterDefaults.(p.fieldname) = p.default;
        x.parameterDefaultValues{i} = p.default;
        x.parameterNames{i} = p.fieldname;
      end
      x.defaultParam = nnetParam(x.mfunction,x.parameterDefaults,x.parameters);
    end
    
  end
    
end

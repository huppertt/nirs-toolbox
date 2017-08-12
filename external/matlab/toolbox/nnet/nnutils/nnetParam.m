classdef nnetParam

% Copyright 2010-2012 The MathWorks, Inc.
  
  properties
    fcn = '';
    values = struct;
    parameters = []; % TODO - Make this lazy evaluation
  end
  
  methods
    
    function x = nnetParam(fcn,values,parameters)
      if isdeployed
        if nargin < 1, fcn = ''; end
        if nargin < 2, values = struct; end
        if nargin < 3, parameters = []; end
        x.fcn = fcn;
        x.values = values;
        x.parameters = parameters;
        return
      end
      if nargin < 1
        fcn = '';
      elseif ~ischar(fcn)
        error(message('nnet:FcnName:NotString'));
      end
      if nargin < 2
        if isempty(fcn)
          values = struct;
        else
          values = feval(fcn,'parameterDefaults');
        end
      elseif isempty(values)
        values = struct;
      elseif isa(values,'nnetParam')
        values = struct(values);
      elseif ~isstruct(values)
        error(message('nnet:nnetParam:NotStruct'));
      end
      if nargin < 3
        if isempty(fcn)
          parameters = [];
        elseif exist(nnpath.fcn2file([fcn '.parameterInfo']),'file')
          parameters = feval([fcn '.parameterInfo']);
        else
          parameters = feval(fcn,'parameters');
        end
      elseif ~isempty(parameters) && ~isa(parameters,'nnetParamInfo')
        error(message('nnet:nnetParam:NotParamInfo'));
      end
      x.fcn = fcn;
      x.values = values;
      x.parameters = parameters;
    end
    
    function disp(x)
      if ~isstruct(x.values)
        error(message('nnet:nnetParam:BadInit'))
      end
      isLoose = strcmp(get(0,'FormatSpacing'),'loose');
      str = {};
      if (isLoose), str{end+1} = ' '; end
      if isempty(x.fcn)
        str{end+1} = '    No Neural Function Parameters';
      elseif isempty(x.parameters)
        str{end+1} = ['    No Function Parameters for ' nnlink.fcn2strlink(x.fcn)];
      else
        str{end+1} = ['    Function Parameters for ' nnlink.fcn2strlink(x.fcn)];
        if (isLoose), str{end+1} = ' '; end
        fields = fieldnames(x.values);
        maxLen1 = 0;
        maxLen2 = 0;
        maxLen3 = 0;
        for i=1:length(fields)
          fi = fields{i};
          ti = x.parameters(i).title;
          maxLen1 = max(maxLen1,length(fi));
          maxLen2 = max(maxLen2,length(ti));
          maxLen3 = max(maxLen3,length(fi) + length(ti));
        end
        for i=1:length(fields)
          fi = fields{i};
          ti = x.parameters(i).title;
          s3 = nnstring.spaces(maxLen3-length(fi)-length(ti));
          stri = ['    ' ti s3 ' ' nnlink.prop2link2(fi) ': '];
          if nnstring.ends(x.parameters(i).type,'_fcn')
            str{end+1} = [stri nnlink.fcn2strlink(x.values.(fi))];
          else
            str{end+1} = [stri nnstring.fieldvalue2str(x.values.(fi))];
          end
        end
      end
      if (isLoose), str{end+1} = ' '; end
      str = nnlink.filterLinks(str);
      for i=1:length(str), disp(str{i}), end
    end

    
    function f = mfunction(x)
      f = x.fcn;
    end
    
    function x = subsref(x,s) % TODO - do I need this?
      x = subsref(x.values,s);
    end
    
    function x = subsasgn(x,s,v)
      x.values = subsasgn(x.values,s,v);
    end
    
    function f = fieldnames(x)
      f = fieldnames(x.values);
    end
    
    function b = isfield(x,fn)
      b = isfield(x.values,fn);
    end
    
    function s = struct(x)
      s = x.values;
    end
    
    function flag = isempty(x)
      flag = isempty(x.parameters);
    end
    
  end
end

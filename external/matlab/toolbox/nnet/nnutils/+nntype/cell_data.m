function [out1,out2] = cell_data(in1,in2,in3)
%NN_CELL_DATA Neural network data in cell array form.

% Copyright 2010 The MathWorks, Inc.

%% =======================================================
%  BOILERPLATE_START
%  This code is the same for all Type Functions.

  persistent INFO;
  if isempty(INFO), INFO = get_info; end
  if nargin < 1, error(message('nnet:Args:NotEnough')); end
  if ischar(in1)
    switch (in1)
      
      case 'info'
        % this('info')
        out1 = INFO;
        
      case 'isa'
        % this('isa',value)
        out1 = isempty(type_check(in2));
        
      case {'check','assert','test'}
        % [*err] = this('check',value,*name)
        nnassert.minargs(nargin,2);
        if nargout == 0
          err = type_check(in2);
        else
          try
            err = type_check(in2);
          catch me
            out1 = me.message;
            return;
          end
        end
        if isempty(err)
          if nargout>0,out1=''; end
          return;
        end
        if nargin>2, err = nnerr.value(err,in3); end
        if nargout==0, err = nnerr.value(err,'Value'); end
        if nargout > 0
          out1 = err;
        else
          throwAsCaller(MException(nnerr.tag('Type',2),err));
        end
        
      case 'format'
        % [x,*err] = this('format',x,*name)
        err = type_check(in2);
        if isempty(err)
          out1 = strict_format(in2);
          if nargout>1, out2=''; end
          return
        end
        out1 = in2;
        if nargin>2, err = nnerr.value(err,in3); end
        if nargout < 2, err = nnerr.value(err,'Value'); end
        if nargout>1
          out2 = err;
        else
          throwAsCaller(MException(nnerr.tag('Type',2),err));
        end
        
      case 'check_param'
        out1 = '';
      otherwise,
        try
          out1 = eval(['INFO.' in1]);
        catch me, nnerr.throw(['Unrecognized first argument: ''' in1 ''''])
        end
    end
  else
    error(message('nnet:Args:Unrec1'))
  end
end

%  BOILERPLATE_END
%% =======================================================

function info = get_info
  info = nnfcnType(mfilename,'NN Cell Data',7.0);
end

function err = type_check(x)
  if ~iscell(x)
    err = 'VALUE is not a cell array.';
  elseif isempty(x)
    err = '';
  elseif ndims(x) ~= 2
    err = 'VALUE is not a 2-dimensional cell array.';
  else
    m = x{1,1};
    
    if ~(isnumeric(m) || islogical(m))
      err = 'VALUE is not numeric or logical.';
    elseif ~isreal(m)
      err = 'VALUE is complex.';
    elseif ndims(m) ~= 2
      err = 'VALUE is not two-dimensional.';
    else
      err = '';
    end
    if ~isempty(err)
      err = nnerr.value(err,'VALUE{1,1}');
      return
    end
    
    [rows,cols] = size(m);
    [cell_rows,cell_cols] = size(x);
    for i=1:cell_rows
      m = x{i,1};
      rows = size(m,1);
      isEmpty = isempty(m);
      for j=1:cell_cols
        m = x{i,j};
        
        if ~(isnumeric(m) || islogical(m))
          err = 'VALUE is not numeric or logical.';
        elseif ~isreal(m)
          err = 'VALUE is complex.';
        elseif ndims(m) ~= 2
          err = 'VALUE is not two-dimensional.';
        else
          err = '';
        end
        if ~isempty(err)
          err = nnerr.value(err,['VALUE{', num2str(i) ',' num2str(j) '}']);
          return;
        end
        
        if isEmpty
          if ~isempty(m)
            err = ['VALUE{' num2str(i) ',1} is empty, but ' ...
              'VALUE{' num2str(i) ',' num2str(j) '} is not.'];
            return;
          end
        else
          [r,c] = size(m);
          if (rows ~= r)
            err = ['VALUE{' num2str(i) ',' num2str(j) '} and VALUE{' ...
              num2str(i) ',1} have different numbers of rows.'];
            return;
          end
          if (cols ~= c)
            err = ['VALUE{' num2str(i) ',' num2str(j) '} and ' ...
              'VALUE{1,1} have different numbers of columns.'];
           return;
          end
        end
      end
    end
    err = '';
  end
end

function x = strict_format(x)
  [S,TS] = size(x);
  Q = 0;
  for i=1:numel(x)
    xi = x{i};
    if ~isempty(xi)
      Q = size(xi,2);
      break
    end
  end
  N = zeros(S,1);
  for i=1:S
    for ts=1:TS
      xits = x{i,ts};
      if ~isempty(xits)
        N(i) = size(xits,1);
        break;
      end
    end
  end
  for i=1:S
    for ts=1:TS
      xits = x{i,ts};
      if isempty(xits)
        x{i,ts} = zeros(N(i),Q);
      elseif ~isa(xits,'double')
        x{i,ts} = double(xits);
      end
    end
  end
end


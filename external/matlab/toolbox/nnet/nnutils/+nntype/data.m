function [out1,out2] = data(in1,in2,in3)
%NN_DATA Neural network data in matrix or cell array form.

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
  info = nnfcnType(mfilename,'NN Data',7.0);
end

function err = type_check(x)
  if iscell(x)
    err = nntype.cell_data('check',x);
  elseif isnumeric(x) || islogical(x)
    err = nntype.matrix_data('check',x);
  else
    % TODO - More detailed response
    err = 'VALUE is not a matrix or cell array.';
  end
end

function x = strict_format(x)
  if islogical(x) || isnumeric(x)
    x = {x};
  else
    [S,TS] = size(x);
    Q = 0;
    Qz = 0;
    for i=1:numel(x)
      xi = x{i};
      if all(size(xi) == 0)
        Q = size(xi,2);
        break;
      elseif size(xi,2) ~= 0
        Qz = size(xi,2);
      end
      if (Q==0), Q = Qz; end
    end
    N = zeros(S,1);
    for i=1:S
      for ts=1:TS
        xits = x{i,ts};
        if all(size(xits) == 0)
          N(i) = size(xits,1);
          break;
        end
      end
    end
    for i=1:S
      for ts=1:TS
        xits = x{i,ts};
        if all(size(xits) == 0)
          x{i,ts} = zeros(N(i),Q);
        elseif ~isa(xits,'double')
          x{i,ts} = double(xits);
        end
      end
    end
  end
end

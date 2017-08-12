function [errorcode,varargout] = distchck(nparms,varargin)
%DISTCHCK Checks the argument list for the probability functions.

%   Copyright 1993-2014 The MathWorks, Inc. 


errorcode = 0;
varargout = varargin;

if nparms == 1
    return;
end

% Get size of each input, check for scalars, copy to output
scalar = cellfun( @isscalar, varargin );

% Done if all inputs are scalars.  Otherwise fetch their common size.
if (all(scalar)), return; end

n = nparms;

sz = cellfun( @size, varargin, 'UniformOutput', false );
t = sz(~scalar);
size1 = t{1};

% Scalars receive this size.  Other arrays must have the proper size.
for j=1:n
   sizej = sz{j};
   if (scalar(j))
      vj = varargin{j};
      if isnumeric(vj)
         t = zeros(size1,'like',vj);
      else
         t = zeros(size1);
      end
      t(:) = vj;
      varargout{j} = t;
   elseif (~isequal(sizej,size1))
      errorcode = 1;
      return;
   end
end

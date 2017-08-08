function varargout = sigprivate(varargin)
%SIGPRIVATE This function allows access to the functions in the private directory.
%        SIGPRIVATE('FOO',ARG1,ARG2,...) is the same as
%        FOO(ARG1,ARG2,...).  
%
%     Copyright 1988-2002 The MathWorks, Inc.

if (nargout == 0)
  feval(varargin{:});
else
  [varargout{1:nargout}] = feval(varargin{:});
end

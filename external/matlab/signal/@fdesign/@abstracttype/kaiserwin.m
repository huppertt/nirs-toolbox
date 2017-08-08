function varargout = kaiserwin(this, varargin)
%KAISERWIN   Design a filter using a kaiser window.
%   KAISERWIN(D) Design a filter using a kaiser window and the
%   specifications in the object D.

%   Author(s): J. Schickler
%   Copyright 1999-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'kaiserwin', varargin{:});
catch e
    throw(e);
end

% [EOF]

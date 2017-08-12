function tag = tag(tag,n)
%NNTAG Updates error/warning tag with NNET and the calling function name.
%
%  NNTAG(TAG) returns 'nnet:FCN:TAG' where FCN is the name
%  of the function which called NNTAG.
%
%  NNTAG(TAG,1) returns the same value.
%
%  NNTAG(TAG,N), where N>1 returns the new tag with FCN being the name
%  of the function N steps up the calling stack.
%
%  If the calling stack is not as long as implied by N then TAG is
%  returned as 'nnet:workspace:TAG';

% Copyright 2010 The MathWorks, Inc.

% Default stack level
if nargin < 2, n = 1; end

% Tag Abbreviations
if strcmpi(tag,'args')
  tag = 'Arguments';
end

% Update Tag
fcn = nnerr.caller(n);
if ~isempty(fcn)
  tag = ['nnet:' fcn ':' tag];
else
  tag = ['nnet:workspace:' tag];
end

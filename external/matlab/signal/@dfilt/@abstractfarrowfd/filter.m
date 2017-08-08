function y = filter(this,x,dim)
%FILTER   

%   Copyright 2005-2012 The MathWorks, Inc.

narginchk(1,3);

if nargin<2, x = []; end
if nargin<3, dim=[]; end

if isempty(x), 
  y = x;
  return; 
end

if ~isnumeric(x) && ~ischar(x) && ~islogical(x)
 error(message('signal:dfilt:abstractfilter:super_filter:InvalidInputType'));
end

s = size(x);
[x,perm,nshifts] = shiftdata(x,dim);
s_shift = size(x); % New size
x = reshape(x,size(x,1),[]); % Force into 2-D

% At this point, x is a 2-D matrix and we always filter along the columns
[Mx,Nx] = size(x);
if log2(Mx*Nx)>31, 
    error(message('signal:dfilt:abstractfarrowfd:filter:InvalidInput'));
end

nchannels = this.nchannels;

if ~this.PersistentMemory,
    % Reset the filter
    reset(this);
else
	if ~isempty(nchannels) && Nx ~= nchannels
		error(message('signal:dfilt:abstractfarrowfd:filter:InvalidDimensions'));
	end
end

% Set number of channels
this.nchannels = Nx;

zi = this.HiddenStates;
% Expand the states for the multichannel case
zi = ziexpand(this,x,zi);

[y,zf] = secfilter(this,x,this.privfracdelay,zi);
this.NumSamplesProcessed = this.NumSamplesProcessed+Mx*Nx;

this.HiddenStates = zf;

if isempty(dim),
    dim = find(s>1,1,'first');
    if isempty(dim), dim = 1; end
end
ly = size(y,1);
s(dim) = ly;
s_shift(1) = ly;
y = reshape(y,s_shift); % Back to N-D array
y = unshiftdata(y,perm,nshifts);

y = reshape(y,s);

% [EOF]

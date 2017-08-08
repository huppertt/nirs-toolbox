function d = defaultmethod(this)
%DEFAULTMETHOD

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

d = designmethods(this, 'fir');
if isempty(d)

    % Only add IIR designs if there are no FIR designs.
    d = designmethods(this, 'iir');
end

if length(d) == 1

    % If there is only 1 available method, use it.
    d = d{1};
elseif any(strcmpi(d, getdefaultmethod(this)))

    % Always use the object specific default first.
    d = getdefaultmethod(this);
elseif any(strcmpi(d, 'equiripple'))

    % Always use equiripple first after the default.
    d = 'equiripple';
elseif any(strcmpi(d, 'ellip'))

    % Always use ELLIPTIC is we have only IIR filters.
    d = 'ellip';
else
    % Otherwise just use the first available one.
    d = d{1};
end


% [EOF]

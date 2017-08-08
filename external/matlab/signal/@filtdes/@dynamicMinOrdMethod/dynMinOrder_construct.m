function dynMinOrder_construct(d, type)
%DYNMINORDER_CONSTRUCT  'Real' constructor for the dynamic min order design method.


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

if nargin < 2,
    type = 'SignalFdatoolFilterOrderMode';
end

% Create a dynamic property to hold the order Mode
% Use a dynamic property because minimum order is not
% always supported (it depends on the filterType in remez for example)
p = schema.prop(d,'orderMode', type);
set(d,'orderMode','minimum');

% Install a listener to the orderMode property
l = handle.listener(d,p,'PropertyPostSet',@filterType_listener);
    
% Store listener
set(l, 'callbacktarget', d); % Allow for methods as callbacks
set(d,'ordModeListener',l);

% Call super's constructor, do this after adding the orderMode prop
singleOrder_construct(d);


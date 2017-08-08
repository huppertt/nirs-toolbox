function designMethodwFs_construct(d,varargin)
%DESIGNMETHODWFS_CONSTRUCT Real constructor for the design method with Fs object.
%
%   Inputs:
%       freqUnits - Frequency units
%       Fs    - Sampling frequency 


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Create a dynamic property to hold the sampling frequency
p = schema.prop(d,'Fs','udouble');
set(d,'Fs',48000);

% Call the super constructor
designMethod_construct(d,varargin{:});

if nargin > 1, set(h,'freqUnits',freqUnits); end
if nargin > 2, set(h,'Fs',Fs); end
    
% Install a listener to the frequency units
l = handle.listener(d,findprop(d,'freqUnits'),'PropertyPreSet',@freqUnits_listener);
    
% Store listener
set(l, 'callbacktarget', d); % Allow for methods as callbacks
set(d,'freqUnitsListener',l);

% Fire the listener manually the first time
freqUnits_listener(d);









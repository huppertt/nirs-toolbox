function propname = determine_dynamicprop(d,freqspec,freqspecOpts)

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


switch freqspec
case freqspecOpts{1}, %'cutoff'
    propname = 'Fc';
case freqspecOpts{2}, %'passedge'
    propname = 'Fpass';
case freqspecOpts{3}, %'stopedge'
    propname = 'Fstop';
end
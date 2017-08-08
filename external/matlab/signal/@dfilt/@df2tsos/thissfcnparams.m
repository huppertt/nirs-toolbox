function varargout = thissfcnparams(Hd)
%THISSFCNPARAMS Returns the parameters for SDSPFILTER

% Author(s): J. Schickler
% Copyright 1988-2002 The MathWorks, Inc.

varargout = {3, 8, getsosstrs(Hd), '', '0'};

% ----------------------------------------------
function sosstrs = getsosstrs(Hd)

sos    = get(Hd, 'sosMatrix');
gain   = get(Hd, 'ScaleValues');
ns     = nsections(Hd);
for indx = 1:min([length(gain) ns]),
    sos(indx, 1:3) = sos(indx, 1:3)*gain(indx);
end

% If there is an output scalevalue do something with it.
if length(gain) > ns,
    sos(end, 1:3) = sos(end, 1:3)*gain(end);
end

sos    = dspblkbiquad2('init', sos);
[r, c] = size(sos);

for indx = 1:r
    
    % Convert each section to a string
    sosstrs{indx} = sprintf('%.25g, ', sos(indx,:));
    
    % Add a semi-colon at the end of each section
    sosstrs{indx}(end-1:end) = ['; '];
end
sosstrs{r}(end-1:end) = [];

% Convert the sections into one string vector
sosstrs = [sosstrs{:}];

% [EOF]

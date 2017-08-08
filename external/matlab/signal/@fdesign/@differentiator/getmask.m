function [F, A] = getmask(this, fcns, ~, specs)
%GETMASK Get the mask.

%   Copyright 2005-2011 The MathWorks, Inc.

w = warning('off'); %#ok<WNOFF>

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
end

apass = specs.Apass;
fpass = specs.Fpass;
astop = specs.Astop;

if isnan(fpass),
    fpass = 1;
end

units = lower(fcns.getunits());
% Get low stopband points
stoplow = fcns.findbottom(-astop);

% Calculate the differentiator frequency vector and amplitude
if strcmpi(units,'dB')
    dc = 10^(stoplow/20);
else
    dc = 0;
end
diffF = linspace(dc, fpass, 100); 
diffA = ones(2,1)*20*log10(diffF*pi); % dB
if ~isnan(apass),
    diffA = diffA+[apass/2;-apass/2]*ones(1,100);
end

offset = sum(diffA(:,end))/2;
stophigh = -astop+offset;

% Convert the differentiator amplitude to the correct units.
if ~strcmpi(units,'dB')
    target = units;
    if strcmpi(units,'zerophase'), target = 'linear'; end
    stoplow = convertmagunits(-stoplow, 'db', target, 'stop');
    stophigh = convertmagunits(-stophigh, 'db', target, 'stop');
    % If we are not given an Astop, we want to draw the Fstop line
    % to 0 when not in dB.
    if isnan(astop),
        stoplow = 0;
        stophigh = 0;
    end
end
stophigh = stophigh*ones(2);
if isnan(stophigh(1,1)),
    stophigh(1,1) = stoplow;
end

switch units
    case 'squared'
        diffA = (10.^(diffA/20)).^2;
    case {'linear', 'zerophase'}
        diffA = 10.^(diffA/20);
        if strcmpi(units, 'zerophase'),
            stoplow = -stoplow;
            stophigh(2,:) = -stophigh(2,:);
        end
end

% Frequency vector 
F = [diffF NaN specs.Fstop specs.Fstop NaN diffF specs.Fpass 1 NaN specs.Fstop 1]*fcns.getfs()/2;

% Construct an amplitude vector.
A = [diffA(1,:) NaN diffA(2,end) stoplow NaN diffA(2,:) stophigh(1,:) NaN stophigh(2,:)];

warning(w);

% [EOF]

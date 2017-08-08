function str = tostring(this)
%TOSTRING

%   Copyright 2005-2009 The MathWorks, Inc.

costfn = fieldnames(this);

for nc = 1:length(costfn),
    value{nc} = num2str(getfield(this,costfn{nc}));
end

param = cell(length(costfn), 1);

for k = 1:length(costfn),
    switch lower(costfn{k}),
        case 'nmult',
            param{k} = getString(message('signal:dfilt:info:NumberofMultipliers'));
        case 'nadd',
            param{k} = getString(message('signal:dfilt:info:NumberofAdders'));
        case 'nstates',
            param{k} = getString(message('signal:dfilt:info:NumberofStates'));
        case 'multperinputsample'
            param{k} = getString(message('signal:dfilt:info:MultiplicationsPerInputSample'));
        case 'addperinputsample'
            param{k} = getString(message('signal:dfilt:info:AdditionsPerInputSample'));
        otherwise
            param{k} = getTranslatedString('signal:dfilt:info',costfn{k});
    end
end
str = [strvcat(param) repmat(' : ', length(param), 1) strvcat(value)];

% [EOF]

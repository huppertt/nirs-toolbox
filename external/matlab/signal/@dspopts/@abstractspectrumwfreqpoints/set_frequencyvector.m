function FrequencyVector = set_frequencyvector(this, FrequencyVector) %#ok
%SET_FREQUENCYVECTOR   PreSet function for the 'FrequencyVector' property.

%   Author(s): W. Syed
%   Copyright 1988-2006 The MathWorks, Inc.

% Welch uses segment length instead of input length.
% auto = max(256,inputlength)

validStrs = {'Auto'};
if ~isnumeric(FrequencyVector),
    
    idx = [];
    for k=1:length(validStrs),
        if regexp(lower(validStrs{k}),['^',lower(FrequencyVector)],'once');
            idx=k;
        end
    end
    
    if isempty(idx),
        error(message('signal:dspopts:abstractspectrumwfreqpoints:set_frequencyvector:invalidFrequencyVectorValue', 'FrequencyVector', 'FrequencyVector', validStrs{ 1 }));
    else
        % Use full string with correct capitalization.
        FrequencyVector = validStrs{idx};
    end
end

% [EOF]

function blk = changeDelayValue(blk,offset,initcond)
%CHANGEDELAYVALUE Offsets delay block's delay value

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.

% get delay latency
delay_str = blk.mainParam;
t = regexpi(delay_str,',');
latency = delay_str(1:t-1);

% create initial condition in the form of column vector (e.g. '[1;2;3]')
init_str = '0';
currentIC = delay_str(t+1:end);
if ~isempty(currentIC)&&~isempty(initcond)
    % rearrange the order of the initial condition. The order of the
    % values in 'initcond' follows the coefficient order of the filter. We
    % need to flip the order so that the last initial condition comes out
    % first from the delay.
    initcond = fliplr(initcond);
    init_str = '[';   
    init_str = sprintf('%s%s',init_str,currentIC);
    for k=1:length(initcond)
        init_str = sprintf('%s;%s',init_str,initcond{k});
    end
    % Need transposition as the delay block consider one row as one
    % channel, which is opposite from the way the filter object stores its
    % state information (each column represents a channel). [g557750]
    currentIC = sprintf('%s].''',init_str); 
end
    
% new delay latency
nvalue = num2str(str2num(latency)+offset);

% restore delay string (letency,initialcondition)
delay_param = [nvalue ',' currentIC];
blk.mainParam = delay_param;
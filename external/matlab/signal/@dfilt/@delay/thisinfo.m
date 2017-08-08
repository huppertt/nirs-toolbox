function [p, v] = thisinfo(this)
%THISINFO   

%   Author(s): M. Chugh
%   Copyright 2005 The MathWorks, Inc.

p = {getString(message('signal:dfilt:info:FilterStructure')), ...
     getString(message('signal:dfilt:info:Latency'))};
v = {get(this, 'FilterStructure'), num2str(get(this,'Latency'))};

if isquantized(this)
    p{end+1} = getString(message('signal:dfilt:info:Arithmetic'));
    v{end+1} = get(this, 'Arithmetic');
        
    [f, c] = info(this.filterquantizer);

    if ~isempty(f)
        p = [p f(:)'];
        v = [v c(:)'];
    end
end

% [EOF]

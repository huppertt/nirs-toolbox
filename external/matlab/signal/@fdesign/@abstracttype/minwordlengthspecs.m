function [Fpass, Fstop, Apass, Astop] = minwordlengthspecs(this,h)
%MINWORDLENGTHSPECS Determine the specs of the minimum wordlength filter

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

% Fpass, Fstop
md = measure(reffilter(h)); % Make sure to use reffilter in case h has been quantized
Fpass = md.Fpass;
Fstop = md.Fstop;

% Astop
if isprop(this,'Astop'),
    if isprop(md,'Astop'),
        % Use as reference in case original (double-precision)
        % design doesn't meet specs.
        Astop = min(md.Astop,this.Astop);
    else
        Astop = this.Astop;
    end
else
    Astop = md.Astop;
end

% Apass
Apass = minwordlengthApass(this,md,Astop);



% [EOF]

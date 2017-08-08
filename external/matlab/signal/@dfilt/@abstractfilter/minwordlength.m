function minwordlength(h,args)
%MINWORDLENGTH Minimum-wordlength design.
%   This method operates in-place so that it can be called from set_arith
%   method. Noise shaping is not used in this method!

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

% Initial estimate of wordlength
wl = ceil(args.Astop/5);

% Quantize filter with initial guess for wordlength
h.privArithmetic = 'fixed';
%h.Arithmetic = 'fixed';
h.CoeffWordLength = wl;

%--------------------------------------------------------------------------
% Find the minimum word length filter that meets the specs
%--------------------------------------------------------------------------
done  = false;
count = 0;
try_smaller = true; % Try smaller wordlengths

% Increase word length until spec is met
while ~done,
    count = count + 1;
    done = isspecmet(h,args.fdesignObj,args);    
    
    if ~done,
        try_smaller = false;
        h.CoeffWordLength = h.CoeffWordLength + 1;
    end
    
    if count == 100,
        error(message('signal:dfilt:abstractfilter:minwordlength:DidNotConverge'));
    end
end

% Decrease word length until spec is no longer met
count = 0;
while try_smaller,
    count = count + 1;
    h.CoeffWordLength = h.CoeffWordLength - 1;
    pass = isspecmet(h,args.fdesignObj,args);
    if ~pass,
        try_smaller = false;
        % Re-increment wordlength
        h.CoeffWordLength = h.CoeffWordLength + 1;
    end    
end

% [EOF]

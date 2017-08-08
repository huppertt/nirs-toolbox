function Hd = thisdesign(d)
%THISDESIGN Design a filter with GREMEZ 

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[F, A, W, args] = getarguments(d.ResponseTypeSpecs, d);

incVals = diff(F);
if any(incVals <= 0)
    error(message('signal:sigtools:fmethod:FrequencySpecMustHaveIncreasingOrder'))
end

if isspecify(d),
    sp  = convertspecialprops(d);
    
    % If there are any special frequency points, insert them before the
    % weights.  Otherwise just add the extra args after the weights.
    if isempty(sp), args = {W, args{:}};
    else,           args = {sp, W, args{:}}; end
    
    % If there are any constraints (IAEs or CEMs) add them to the end.
    con = getconstraints(d);
    if ~isempty(con), args = {args{:}, con}; end
    
else
    args = {W, args{:}};
end

order = convertorder(d);

args = {args{:}, {get(d, 'DensityFactor')}};

phase = get(d, 'Phase');
if ~strcmpi(phase, 'linear'),
    args = {args{:}, sprintf('%sphase', lower(phase(1:3)))};
end

if isdynpropenab(d, 'FIRType'),
    firtype = get(d, 'FIRType');
    if ~strcmpi(firtype, 'unspecified'), args = {args{:}, lower(firtype)}; end
end

try
    
    b = feval(designfunction(d), order, F, A, args{:});

catch ME
    
    % In min order we cannot error about filter order being too large
    if strcmp(d.orderMode,'minimum') && strcmp(ME.identifier,'dsp:firminphase:InvalidDimensions')
        error(message('signal:sigtools:fmethod:TryRelaxingSomeSpecifications'))        
    end
    
end


% Construct object
Hd = dfilt.dffir(b);






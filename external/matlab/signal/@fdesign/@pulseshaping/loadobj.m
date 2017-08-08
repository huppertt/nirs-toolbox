function this = loadobj(s)
%LOADOBJ   Load this object.

%   Copyright 2008 The MathWorks, Inc.

this = feval(s.class);

this.PulseShape    = s.PulseShape;

if ~isfield(s, 'version')
    % This was saved with 8b or before
    
    % First remove 'PulseShape" from specs, if exists
    AllSpecs = s.AllSpecs;
    for p=1:length(AllSpecs)
        if isfield(AllSpecs{p}, 'PulseShape')
            AllSpecs{p} = rmfield(AllSpecs{p}, 'PulseShape');
        end
    end
    s.AllSpecs = AllSpecs;
    
    loadobj(this.PulseShapeObj, s);
    
    % If the old spec was square root raised cosine, that spec did not exist
    % independently in 8b, instead was a part of fspecs.psrcosnsym.  We need to
    % update the properties manually, since updatecurrentspecs did not recognize
    % fspecs.psrcosnsym in AllSpecs and created a new fspecs.pssqrtrcosnsym.
    cSpecCon = getconstructor(this.PulseShapeObj);
    idx = strmatch(cSpecCon, {'fspecs.pssqrtrcosnsym', 'fspecs.pssqrtrcosord'});
    if ~isempty(idx)
        if idx == 1
            for p=1:length(AllSpecs)
                if strcmp(AllSpecs{p}.class, 'fspecs.psrcosnsym')
                    spec = AllSpecs{p};
                end
            end
            this.NumberOfSymbols = spec.NumberOfSymbols;
        else
            for p=1:length(AllSpecs)
                if strcmp(AllSpecs{p}.class, 'fspecs.psrcosord')
                    spec = AllSpecs{p};
                end
            end
            this.FilterOrder = spec.FilterOrder;
        end
        this.NormalizedFrequency = spec.NormalizedFrequency;
        if this.NormalizedFrequency == false
            this.Fs = spec.Fs;
        end
        this.SamplesPerSymbol = spec.SamplesPerSymbol;
        this.RolloffFactor = spec.RolloffFactor;
    end
else
    % This was saved with 9a or after
    loadobj(this.PulseShapeObj, s.PulseShapeObj);
end

% [EOF]

function Href = optimizestopbandfir(this,Href,wl,varargin) %#ok<INUSL>
%OPTIMIZESTOPBANDFIR Optimize stopband.
%   This should be a private method

%   Author(s): R. Losada
%   Copyright 2009 The MathWorks, Inc.

hf = getfdesign(Href);
hm = getfmethod(Href);

args = sboptiminpututparse(Href,varargin{:});


% make a copy of input obj so we don't change it
Href.privArithmetic = 'fixed';
Href.CoeffWordLength = wl;
m = measure(Href);

bestk = 0;
bestAstop = m.Astop;
for k = 1:args.NTrials,
    Hns(k) = noiseshape(hf,Href,wl,args); %#ok<AGROW>
    m(k) = measure(Hns(k),hf);
    if m(k).Astop > bestAstop,
        bestAstop = m(k).Astop;
        bestk = k;
    end
end

if bestk > 0,
    % We found a better filter
    Href = Hns(bestk);
    Href.privArithmetic = 'fixed';
    Href.CoeffWordLength = wl;
        
    % Set fdesign, fmethod so that FVTool shows specs
    setfdesign(Href, hf);
    setfmethod(Href, hm);
end

%--------------------------------------------------------------------------
function args = sboptiminpututparse(this,varargin)

% Test if response type is supported
s = iscoeffwloptimizable(this);

p = inputParser; 
p.addParamValue('NTrials',1,@(x)(isnumeric(x)&&x>=1&&x==round(x)));
p.addParamValue('noiseShapeNorm',inf,@(x)(isnumeric(x)&&(x>0)));
p.parse(varargin{:});
args = p.Results;

% Copy fieds of s into args
f = fieldnames(s);
for i=1:length(f),
    args.(f{i}) = s.(f{i});
end

% Use reference filter in case input filter has been quantized
hmethod        = getfmethod(this);
method         = hmethod.DesignAlgorithm;
args.dm        = method;
args.dopts     = designopts(hmethod);





% [EOF]

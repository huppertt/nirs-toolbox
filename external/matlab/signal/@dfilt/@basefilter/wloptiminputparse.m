function args = wloptiminputparse(this,varargin)
%WLOPTIMINPUTPARSE Parse inputs for wordlength optimization functions.

% This should be a private method

%   Copyright 2009 The MathWorks, Inc.

% Test if response type is supported
s = iscoeffwloptimizable(this);

p = inputParser; 
p.addParamValue('NTrials',1,@(x)(isnumeric(x)&&x>=1&&x==round(x)));
p.addParamValue('Apasstol',1e-4,@(x)(isnumeric(x)&&x>=0));
p.addParamValue('Astoptol',1e-2,@(x)(isnumeric(x)&&x>=0));
p.addParamValue('noiseShaping',true,@islogical);
p.addParamValue('noiseShapeNorm',inf,@(x)(isnumeric(x)&&(x>0)));
p.addParamValue('MatchRefFilter',false,@islogical);
p.parse(varargin{:});
args = p.Results;

% Copy fieds of s into args
f = fieldnames(s);
for i=1:length(f),
    args.(f{i}) = s.(f{i});
end

% Add design method and design options
hmethod        = getfmethod(this);
method         = hmethod.DesignAlgorithm;
args.dm        = method;
args.dopts     = designopts(hmethod);



function s = designopts(this, designmethod, sigonlyflag)
%DESIGNOPTS   Return information about the design options.

%   Copyright 2005-2013 The MathWorks, Inc.

narginchk(2,3);

if nargin == 3 
    if strcmp(sigonlyflag,'signalonly')
        s = designopts(this.CurrentSpecs, lower(designmethod),true);   
        if isfield(s,'SystemObject')
            s = rmfield(s,'SystemObject');            
        end
        if isfield(s,'FilterStructure')
            s = rmfield(s,'FilterStructure');            
        end        
        if isfield(s,'DensityFactor')
            s = rmfield(s,'DensityFactor');            
        end 
        if isfield(s,'SOSScaleNorm')
          s = rmfield(s,'SOSScaleNorm');
        end
        if isfield(s,'SOSScaleOpts')
          s = rmfield(s,'SOSScaleOpts');
        end        
    else
        error(message('signal:fdesign:abstracttypewspecs:designopts:invalidSigFlag'))
    end
else
    s = designopts(this.CurrentSpecs, lower(designmethod));
end


% [EOF]

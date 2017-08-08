function constr(this,varargin)
%CONSTR   Constructor for cascade allpass

%   Author(s): R. Losada
%   Copyright 2005-2006 The MathWorks, Inc.

if nargin > 1,    
    for k = 1:length(varargin),
        if ~isreal(varargin{k}),
            error(message('signal:dfilt:abstractcascadeallpass:constr:complexCoeffs'));
        end
        flds{k} = sprintf('Section%d',k);
    end
else
    flds = 'Section1';
    varargin={[]};
end

this.AllpassCoefficients = cell2struct(varargin,flds,2);


% [EOF]

classdef DiscriminantParams < classreg.learning.modelparams.ModelParams
%DiscriminantParams Parameters for discriminant analysis.
%
%   DiscriminantParams properties:
%       DiscrimType       - 'Linear', 'pseudoLinear', 'diagLinear',
%                           'quadratic', 'pseudoQuadratic', or
%                           'diagQuadratic'.
%       Gamma             - Parameter for regularizing the correlation
%                           matrix of predictors.
%       Delta             - Threshold for linear coefficients.
%       FillCoeffs        - Flag for computing discriminant coefficients.
%       SaveMemory        - Flag for deferring computation of the full
%                           covariance matrix. 
    
%   Copyright 2011 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public,Hidden=true)
        AllowedDiscrimTypes = {'linear' 'quadratic' ...
            'diagLinear' 'diagQuadratic' 'pseudoLinear' 'pseudoQuadratic'};
    end

    properties
        DiscrimType = '';
        Gamma = [];
        Delta = [];
        FillCoeffs = [];
        SaveMemory = [];
    end
    
    methods(Access=protected)
        function this = DiscriminantParams(mode,gamma,delta,fillcoeffs,savememory)
            this = this@classreg.learning.modelparams.ModelParams('Discriminant','classification');
            this.DiscrimType = mode;
            this.Gamma = gamma;
            this.Delta = delta;
            this.FillCoeffs = fillcoeffs;
            this.SaveMemory = savememory;
        end
    end

    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin)
            % Decode input args
            args = {'discrimtype' 'gamma' 'delta' 'fillcoeffs' 'savememory'};
            defs = {           ''      []      []           []           []};
            [mode,gamma,delta,fillcoeffs,savememory,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check discriminant type
            allowed = classreg.learning.modelparams.DiscriminantParams.AllowedDiscrimTypes;
            if ~isempty(mode) 
                if ~ischar(mode)
                    error(message('stats:classreg:learning:modelparams:DiscriminantParams:make:DiscrimTypeNotChar'));
                end
                tf = strncmpi(mode,allowed,length(mode));
                if sum(tf)~=1
                    error(message('stats:classreg:learning:modelparams:DiscriminantParams:make:BadDiscrimType', sprintf( ' %s', allowed{ : } )));
                end
                mode = allowed{tf};
            end
            
            % Check gamma
            if ~isempty(gamma) && (~isnumeric(gamma) || ~isscalar(gamma) || gamma<0 || gamma>1)
               error(message('stats:classreg:learning:modelparams:DiscriminantParams:make:BadGamma'));
            end
            
            % Check delta
            if ~isempty(delta) && (~isnumeric(delta) || ~isscalar(delta) || delta<0)
               error(message('stats:classreg:learning:modelparams:DiscriminantParams:make:BadDelta'));
            end
            
            % Check fill coeffs
            if ~isempty(fillcoeffs)
                if ~ischar(fillcoeffs) || (~strcmpi(fillcoeffs,'on') && ~strcmpi(fillcoeffs,'off'))
                    error(message('stats:classreg:learning:modelparams:DiscriminantParams:make:BadFillCoeffs'));
                end
                fillcoeffs = strcmpi(fillcoeffs,'on');
            end
            
            % Check save memory
            if ~isempty(savememory)
                if ~ischar(savememory) || (~strcmpi(savememory,'on') && ~strcmpi(savememory,'off'))
                    error(message('stats:classreg:learning:modelparams:DiscriminantParams:make:BadSaveMemory'));
                end
                savememory = strcmpi(savememory,'on');
            end

            % Make holder
            holder = classreg.learning.modelparams.DiscriminantParams(mode,gamma,delta,fillcoeffs,savememory);
        end
    end

    methods(Hidden)
        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary)
            if isempty(this.DiscrimType)
                this.DiscrimType = 'linear';
            end
            if isempty(this.Gamma)
                this.Gamma = 0;
            end
            if isempty(this.Delta)
                this.Delta = 0;
            end
            if isempty(this.FillCoeffs)
                this.FillCoeffs = true;
            end
            if isempty(this.SaveMemory)
                this.SaveMemory = false;
            end
        end
    end

end

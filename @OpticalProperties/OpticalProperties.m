classdef OpticalProperties
    %OPTICALPROPERTIES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mua;
        mus;
        ri;
        lambda;
    end
    
    properties( Dependent = true )
        kappa;
        c;
    end
    
    methods
        %% Constructor
        function obj = OpticalProperties( varargin )
            if nargin > 0
                obj.mua = varargin{1};
            end
            
            if nargin > 1
                obj.mus = varargin{2};
            end
            
            if nargin > 2
                obj.ri = varargin{3};
            end
            
            if nargin > 3
                obj.lambda = varargin{4};
            end
            
            if nargin > 4
                error( 'Too many arguments.' )
            end
        end
        
        %% Set/Get
        function obj = set.mua( obj, newMua )
            obj.mua = newMua;
            if any(newMua < 0)
                warning( 'Absorption coefficient cannot be negative.' )
            end
        end
        
        function obj = set.mus( obj, newMus )
            obj.mus = newMus;
            if any(newMus < 0)
                warning( 'Scattering coefficient cannot be negative.' )
            end
        end
        
        function obj = set.kappa( obj, newKappa )
            if ~isempty( obj.mua )
                    obj.mus = 1./newKappa/3 - obj.mua;
            end
            
            if any(newKappa < 0)
                warning( 'Diffusion coefficient cannot be negative.' )
            end
        end
        
        function obj = set.ri( obj, newRI )
            if all(newRI >= 1)
                obj.ri = newRI;
            else
                error( 'Refractive index must be greater than or equal to 1.' )
            end
        end
        
        function obj = set.lambda( obj, newLambda )
            if all(newLambda > 0)
                obj.lambda = newLambda;
            else
                error( 'Wavelength must be greater than 0.' )
            end
        end
        
        function kappa = get.kappa( obj )
            kappa = 1 ./ ( 3*obj.mua + 3*obj.mus );
        end
        
        function c = get.c( obj )
           c = 3e11 ./ obj.ri;
        end
        
        %% Methods
        function out = isValid( obj )
            out = length( obj.mua ) == length( obj.lambda );
            out = out && length( obj.mus ) == length( obj.lambda );
            out = out && length( obj.ri ) == length( obj.lambda );
        end
        
    end
end
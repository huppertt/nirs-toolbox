classdef OpticalProp
    %OPTICALPROP This class holds optical properties in terms of absorption
    % and scattering.
    
    properties
        lambda;      	% wavelengths (nm)
        
        mua          	% absorption (mm^-1)
        mus;          	% reduced scattering (mm^-1)
         
        ri;            	% refractive index
    end
    
    properties( Dependent )
        kappa;       	% diffusion coefficient ( 1/(3*mua + 3*mus) )
        v;            	% speed of light in medium
    end
    
    methods
        %% Constructor
        function obj = OpticalProperties( mua, mus, lambda, ri )
            if nargin == 0
                obj.lambda = [690 830];
                obj.mua = [0.0113173914928 0.0133132051072];	
                obj.mus = [1.440357033817907 1.069599610609434];
                obj.ri = 1.45;
            end
            if nargin > 2
                obj.mua = mua;
                obj.mus = mus;
                obj.lambda = lambda;
                obj.ri = 1.45;
            end
            if nargin == 4
                obj.ri = ri;
            end
        end
        
        %% Set Methods
        function obj = set.mua( obj, mua )
            assert( isvector(mua) && all(mua >= 0) )
            obj.mua = mua(:)';
        end
        
        function obj = set.mus( obj, mus )
            assert( isvector(mus) && all(mus >= 0) )
            obj.mus = mus(:)';
        end

        function obj = set.lambda( obj, lambda )
            assert( isvector(lambda) && all(lambda > 0) )
          	obj.lambda = lambda(:)';
        end
        
        function obj = set.ri( obj, ri )
            assert( isscalar(ri) && ri >= 1 );
            obj.ri = ri(:)';
        end
        
        %% Dependent Methods
        function v = get.v( obj )
           v = 3e11 ./ obj.ri;
        end
        
        function kappa = get.kappa( obj )
            kappa = 1./(3*obj.mua + 3*obj.mus);
        end
               
        function obj = set.kappa( obj, kappa )
            assert( isvector(kappa) && kappa >= 0 );
            obj.mus = 1./kappa/3 - obj.mua;
        end
        
        %% Other Methods
        function out = isValid( obj )
            out = length( obj.mua ) == length( obj.lambda ) ...
                && length( obj.mus ) == length( obj.lambda ) ...
                && length( obj.ri ) == length( obj.lambda );
        end
        
    end
end

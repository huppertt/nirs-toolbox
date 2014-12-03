classdef Properties
    %PROPERTIES This class holds optical properties.  Absorption can
    %be parameterized in terms of mua; hbo and hbr; or so2 and hbt.
    %Scattering can be parameterized in terms of mus or a, b such that 
    % mus = a*(lambda/500)^-b
    
    properties
        lambda;                 % wavelengths (nm)
        
        mua                     % absorption (mm^-1)
        mus;                    % reduced scattering (mm^-1)
         
        ri;                     % refractive index
        
        water;                  % water content (0 to 1)
        lipid;
        cytC;
    end
    
    properties( Dependent )
        hbo;                 	% oxy-hemoglobin (uM)
        hbr;                  	% deoxy-hemoglobin (uM)
        
        so2;                 	% oxygen saturation (0 to 1)
        hbt;                 	% total hemoglobin (uM)

        kappa;                  % diffusion coefficient ( 1/(3*mua + 3*mus) )
        
        a;                      % mie scattering magnitude
        b;                  	% mie scaattering exponent
        
        v;                  	% speed of light in medium
    end
    
    methods
        %% Constructor
        function obj = Properties( mua, mus, lambda, ri )
            if nargin == 0
                obj.lambda = [690 830]';
                obj.mua = [0.0113173914928 0.0133132051072]';	
                obj.mus = [1.440357033817907 1.069599610609434]';
                obj.ri = [1.45 1.45]';
                obj.water = 0.7;
                obj.lipid = 0;
                obj.cytC = 0;
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
        
        %% Standard Set Methods
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
            assert( isvector(ri) && all(ri >= 1) );
            obj.ri = ri(:)';
        end
        
        %% Dependent Methods
        function c = get.v( obj )
           c = 3e11 ./ obj.ri;
        end
        
        function kappa = get.kappa( obj )
            kappa = 1./(3*obj.mua + 3*obj.mus);
        end
        
        function obj = set.kappa( obj, kappa )
            assert( isvector(kappa) && kappa >= 0 );
            obj.mus = 1./kappa/3 - obj.mua;
        end
        
        %% Spectral Get Methods
        function a = get.a( obj )
            assert( length(obj.lambda) > 1 )
            
            n = size(obj.lambda,2);
            tmp = pinv( [ones(n,1) -log(obj.lambda'/500)] ) * log( obj.mus )';
            
            a = exp(tmp(1));
        end
        
        function b = get.b( obj )
            assert( length(obj.lambda) > 1 )
            
            n = size(obj.lambda,2);
            tmp = pinv( [ones(n,1) -log(obj.lambda'/500)] ) * log( obj.mus )';
            
            b = tmp(2);
        end

        function hbo = get.hbo( obj )
            assert( length(obj.lambda) > 1 )
            
            ext = nirs2.getSpectra( obj.lambda );
            
            hb = pinv( ext(:,1:2) )  ...
                *  (obj.mua' - ext(:,3)*obj.water + ext(:,4)*obj.lipid + ext(:,5)*obj.cytC);
            
            hbo = roundn( 1e6*hb(1),-12 ); % rounding prevents really small negative numbers
        end
        
        function hbr = get.hbr( obj )
            assert( length(obj.lambda) > 1 )
            
            ext = nirs2.getSpectra( obj.lambda );
            
            hb = pinv( ext(:,1:2) )  ...
                *  (obj.mua' - ext(:,3)*obj.water + ext(:,4)*obj.lipid + ext(:,5)*obj.cytC);
            
            hbr = roundn( 1e6*hb(2),-12 ); % rounding prevents really small negative numbers
        end
        
        function so2 = get.so2( obj )
            assert( length(obj.lambda) > 1 )
            
            so2 = obj.hbo / obj.hbt;
        end
        
        function hbt = get.hbt( obj )
            assert( length(obj.lambda) > 1 )
            
            hbt = obj.hbo + obj.hbr;
        end
        
        %% Spectral Set Methods (only work if >1 wavelength)
        function obj = set.hbo( obj,hbo )
            assert( hbo >= 0 && length(obj.lambda) > 1 )
            
            ext = nirs2.getSpectra( obj.lambda );

            mua = ext(:,1)*hbo*1e-6 + ext(:,2)*obj.hbr*1e-6 + ...
                    ext(:,3)*obj.water + ext(:,4)*obj.lipid + ext(:,5)*obj.cytC;

            obj.mua = mua(:)';
        end
        
        function obj = set.hbr( obj,hbr )
            assert( hbr >= 0 && length(obj.lambda) > 1 )
            
            ext = nirs2.getSpectra( obj.lambda );

            mua = ext(:,1)*obj.hbo*1e-6 + ext(:,2)*hbr*1e-6 + ...
                    ext(:,3)*obj.water + ext(:,4)*obj.lipid + ext(:,5)*obj.cytC;

            obj.mua = mua(:)';
        end
        
        function obj = set.so2( obj,so2 )
            assert( so2 >= 0 && so2 <= 1 && length(obj.lambda) > 1 )
            
            % very important to calculate first, otherwise changing hbo
            % will change hbt
            hbo = so2 * obj.hbt;
            hbr = (1-so2) * obj.hbt;
            
            obj.hbo = hbo;
            obj.hbr = hbr;
        end
        
        function obj = set.hbt( obj,hbt )
            assert( hbt >= 0 && length(obj.lambda) > 1 )
            
            % very important to calculate first, otherwise changing hbo
            % will change hbt
            hbo = obj.so2 * hbt;
            hbr = (1-obj.so2) * hbt;
            
            obj.hbo = hbo;
            obj.hbr = hbr;
        end

        function obj = set.a( obj,a )
            assert( a > 0 && length(obj.lambda) > 1 )
            obj.mus = a * (obj.lambda/500).^-obj.b;
        end
        
        function obj = set.b( obj,b )
            assert( a > 0 && length(obj.lambda) > 1 )
            obj.mus = obj.a * (obj.lambda/500).^-b;
        end
        
        %% Methods
        function out = isValid( obj )
            out = length( obj.mua ) == length( obj.lambda ) ...
                && length( obj.mus ) == length( obj.lambda ) ...
                && length( obj.ri ) == length( obj.lambda );
        end
        
    end
end
classdef SpectralProp
    
    %SPECTRALPROP This class holds optical properties with additional 
    %set/get method for spectral properties.  Absorption can
    %be parameterized in terms of mua; hbo and hbr; or so2 and hbt.
    %Scattering can be parameterized in terms of mus or a, b such that 
    % mus = a*(lambda/500)^-b
    
    properties
        so2;        	% oxygen saturation (0 to 1)
        hbt;          	% total hemoglobin (uM)
        
        water;      	% water content (0 to 1)
        lipid;
        % cytC;
        
        a;          	% mie scattering magnitude
        b;              % mie scaattering exponent 
        
        lambda;         % wavelengths (nm)
        g;              % anisotropy index
        ri;             % refractive index
    end
    
    properties( Dependent )
        mua;                    % absorption (mm^-1)
        mus;                    % scattering (mm^-1)
        musp;                   % reduced scattering (mm^-1)
        hbo;                 	% oxy-hemoglobin (uM)
        hbr;                  	% deoxy-hemoglobin (uM)
   
        kappa;                  % diffusion coefficient ( 1/(3*mua + 3*mus) )
        v;                      % speed of light in medium
    end
    
    methods
        %% Constructor
        function obj = SpectralProp( so2, hbt, lambda, mus, ri,g )
            obj.lambda = [690 830];
            obj.so2 = 0.7;
            obj.hbt = 60;

            obj.a = 2.42;  % mm -1
            obj.b = 1.611;
            obj.g = .89;
            obj.ri = 1.45;
            obj.water = 0.7;
            obj.lipid = 0;
            %obj.cytC = 0;
                
            if nargin > 2
                obj.so2 = so2;
                obj.hbt = hbt;
                obj.lambda = lambda;
                obj.ri = 1.45;
            end
            
            if nargin > 3
                obj.mus = mus;
            end
            
            if nargin > 4
                obj.ri = ri;
            end
            if nargin == 5
                obj.g = g;
            end
        end
        
       %% Set Methods
       function obj = set.lambda( obj, lambda )
          assert( isvector(lambda) && length(lambda) > 0 && all(lambda > 650) && all(lambda < 900) )
          obj.lambda = lambda;
       end
       
       function obj = set.so2( obj, so2 )
           assert( isscalar(so2) && so2 >= 0 && so2 <= 1 )
           obj.so2 = so2;
       end
        function obj = set.g(obj,g)
            assert( isvector(g) && all(g >= 0) )
            obj.g = g(:)';
        end
        
       function obj = set.hbt( obj, hbt )
           assert( isscalar(hbt) && hbt >= 0 )
           obj.hbt = hbt;
       end
       
       function obj = set.water( obj, water )
           assert( isscalar(water) && water >= 0 && water <= 1 )
           obj.water = water;
       end
       
       function obj = set.lipid( obj, lipid )
           assert( isscalar(lipid) && lipid >= 0 )
           obj.lipid = lipid;
       end
       
%        function obj = set.cytC( obj, cytC )
%            assert( isscalar(cytC) && cytC >= 0 )
%            obj.cytC = cytC;
%        end
       
       function obj = set.a( obj, a )
           assert( isscalar(a) && a >= 0 )
           obj.a = a;
       end
       
       function obj = set.b( obj, b )
           assert( isscalar(b) )
           obj.b = b;
       end
       
       function obj = set.ri( obj, ri )
           assert( isscalar(ri) && ri >= 1 )
           obj.ri = ri;
       end
        
        %% Get Methods
        function hbo = get.hbo( obj )
            hbo = obj.so2 * obj.hbt;
        end
        
        function hbr = get.hbr( obj )
           hbr = (1-obj.so2) * obj.hbt;
        end
        
        function mua = get.mua( obj )
            ext = nirs.media.getspectra( obj.lambda );

            mua = ext(:,1)*obj.hbo*1e-6 + ext(:,2)*obj.hbr*1e-6 + ...
                    ext(:,3)*obj.water + ext(:,4)*obj.lipid; % + ext(:,5)*obj.cytC;
                
            mua = mua;
        end
        
        function mus = get.mus( obj )
            mus = (obj.a * (obj.lambda/500).^-obj.b)/(1-obj.g);
        end
        
        function musp = get.musp( obj )
            musp = (obj.a * (obj.lambda/500).^-obj.b);
        end
        
        function v = get.v( obj )
           v = 3e11 ./ obj.ri;
        end
        
        function kappa = get.kappa( obj )
            kappa = 1./(3*obj.mua + 3*obj.musp);
        end
        
        %% Dependent Set Methods
        function obj = set.hbo( obj, hbo )
            hbr = obj.hbr;
            obj.hbt = hbo + hbr;
            obj.so2 = hbo / (hbo + hbr);
        end
        
        function obj = set.hbr( obj, hbr )
            hbo = obj.hbo;
            obj.hbt = hbo + hbr;
            obj.so2 = hbo / (hbo + hbr);
        end
        
        function obj = set.mua( obj, mua )
            assert( isvector(mua) && length(mua) == length(obj.lambda) )
            
            ext = nirs.media.getspectra( obj.lambda );
            hb = pinv( ext(:,1:2) )  ...
                *  (mua(:) - ext(:,3)*obj.water + ext(:,4)*obj.lipid);% + ext(:,5)*obj.cytC);
            
            obj.so2 = hb(1) / sum(hb);
            obj.hbt = sum(hb) * 1e6;
        end
        
        function obj = set.mus( obj, mus )
            assert( isvector(mus) && length(mus) == length(obj.lambda) )
            musp = mus*(1-obj.g);
            n = size(obj.lambda,2);
            p = pinv( [ones(n,1) -log(obj.lambda'/500)] ) * log( musp(:) );
            
            obj.a = exp(p(1));
            obj.b = p(2);
        end
        
        function obj = set.kappa( obj, kappa )
            assert( isvector(kappa) && all(kappa) >= 0 );
            obj.mus = (1./kappa/3 - obj.mua)/(1-obj.g);
        end
        
    end
end

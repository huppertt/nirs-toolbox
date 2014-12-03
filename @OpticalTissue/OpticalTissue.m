classdef OpticalTissue %< nirs.OpticalProperties
    %OPTICALTISSUE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        HbO = 40e-6;
        HbR = 10e-6;
        water = 0.7;
        lipid = 0;
        cytC = 0;

        mus;
        ri = 1.45;
    end
    
    properties( Dependent = true )
        kappa;
        c;
    end

%     methods( Static )
%       	out = getExtinctions( lambda, iSpectrum )
%     end
    
    methods
        %% Constructor
        function obj = OpticalTissue( varargin )
            if nargin > 0
                obj.HbO = varargin{1};
            end
            
            if nargin > 1
                obj.HbR = varargin{2};
            end
            
            if nargin > 2
                obj.mus = varargin{3};
            end
            
            if nargin > 3
                obj.ri = varargin{4};
            end
            
            if nargin > 4
                error( 'Too many arguments.' )
            end
        end
        
        %% Set/Get
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
        
        function kappa = get.kappa( obj )
            kappa = 1 ./ ( 3*obj.mua + 3*obj.mus );
        end
        
        function c = get.c( obj )
           c = 3e11 ./ obj.ri;
        end
        
        %% Methods
        function out = mua( obj, lambda )
            if isvector( lambda )
                if isrow(lambda)
                    lambda = lambda';
                end
                
%                 ext = nirs.getExtinctions( lambda, 1 ) / 10;
                ext = nirs.getSpectra(lambda);

                mua = ext(:,1)*obj.HbO + ext(:,2)*obj.HbR + ...
                    ext(:,3)*obj.water + ext(:,4)*obj.lipid + ext(:,5)*obj.cytC;
                out = mua; % convert log base and convert to mm^-1
            else
                error('Lambda should be vector of wavelenths in nm.')
            end
        end
    end
    
end


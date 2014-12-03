classdef LayeredSpectralModel
    %MCXSPECTRAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        optTissue; 
    end
    
    properties( Dependent = true )
        image;
        probe;
        
        modFreq = 110e6;
        nPhotons = 1e7;
        nRepetitions = 1;
        nTimeGates = 32;
        timeStep = 1/110e6/32;
    end
    
    properties( SetAccess = private )
        nLayers = 1;
    end
    
    properties
        fwdModel;
    end
    
    methods
        %% Constructor 
        function obj = LayeredSpectralModel( varargin )
            obj.fwdModel = nirs.MCXForwardModel();
            
            if nargin > 0
                obj.image = varargin{1};
            end
            
            if nargin > 1
                obj.optTissue = varargin{2};
            end
            
            if nargin > 2
                obj.probe = varargin{3};
            end
            
            if nargin > 3
                obj.modFreq = varargin{4};
            end
            
            if nargin > 4
                error( 'Too many input arguments.' )
            end 
        end
        
        %% Set/Get
        function obj = set.image( obj, newImage )
            if isa( newImage,'nirs.Image' ) && isa( newImage.volume,'uint8' )
                obj.fwdModel.image = newImage;
                obj.nLayers = max( obj.fwdModel.image.volume(:) );
            else
                error('Image should be of Image class with uint8 volume.')
            end
        end
        
        function obj = set.probe( obj, newProbe )
            if isa( newProbe,'nirs.Probe' )
                obj.fwdModel.probe = newProbe;
            else
                error('Probe should be of Probe class.')
            end
        end
        
        function obj = set.optTissue( obj, newTissue )
            if isa( newTissue,'nirs.OpticalTissue' )
                obj.optTissue = newTissue;
            else
                error('optTissue should be of OpticalTissue class.')
            end
        end
        
        function out = get.image( obj )
            out = obj.fwdModel.image;
        end
        
        function out = get.probe( obj )
            out = obj.fwdModel.probe;
        end
        
        %% Methods
        function meas = measurement( obj )
            obj = obj.setProperties();
            meas = obj.fwdModel.measurement();
        end
        
        function [Js,meas] = jacobian( obj )
            obj = obj.setProperties();
            [J, meas] = obj.fwdModel.jacobian();
            
            % convert jacobian
            Js.HbO = J.mua;
            Js.HbR = J.mua;
            Js.water = J.mua;
            Js.lipid = J.mua;
            Js.cytC = J.mua;
            Js.mus = J.kappa;
                        
            for iLambda = 1:length( obj.probe.lambda )
                lst = obj.probe.link(:,3) == iLambda;
                
                for iLayer = 1:obj.nLayers
                    ext = log(10) * nirs.getExtinctions( obj.probe.lambda(iLambda) ) / 10;
                    
                    Js.HbO(lst,iLayer) = Js.HbO(lst,iLayer) * ext(1);
                    Js.HbR(lst,iLayer) = Js.HbR(lst,iLayer) * ext(2);
                    Js.water(lst,iLayer) = Js.water(lst,iLayer) * ext(3);
                    Js.lipd(lst,iLayer) = Js.lipid(lst,iLayer) * ext(4);
                    Js.cytC(lst,iLayer) = Js.cytC(lst,iLayer) * ext(5);
                    
                    Js.mus(lst,iLayer) = -1/3/obj.optTissue(iLayer).mus^2 * J.kappa(lst,iLayer);
                    
                end
            end
        end
    end
    
    methods( Access = private )
        function obj = setProperties( obj )
            lambda = obj.probe.lambda;
            
            optProp(1:obj.nLayers) = nirs.OpticalProperties();
            
            for i = 1:obj.nLayers
                mua = obj.optTissue(i).mua( lambda );
                mus = obj.optTissue(i).mus*ones(size(lambda));
                ri = obj.optTissue(i).ri;
                optProp(i) = nirs.OpticalProperties( mua, mus, ri, lambda );
            end
            
            obj.fwdModel.optProp = optProp;
        end
    end
    
end


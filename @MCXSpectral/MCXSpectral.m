classdef MCXSpectral
    %MCXSPECTRAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        optTissue; 
    end
    
    properties( Dependent = true )
        gpuId;
        image;
        probe;
        
        modFreq = 110e6;
        nPhotons = 1e7;
        nRepetitions = 1;
        nTimeGates = 32;
        timeStep = 1/110e6/32;
    end
    
    properties( SetAccess = private )
        nLayers;
    end
    
    properties( Access = private )
        fwdModel;
    end
    
    methods
        %% Constructor 
        function obj = MCXSpectral( varargin )
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
        function obj = set.gpuId( obj, newID )
            if isnumeric( newID )
                obj.fwdModel.gpuId = newID;
            else
                error('gpuID should be numeric.')
            end
        end
        
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
        
        function obj = set.modFreq( obj, newModFreq )
            if isnumeric( newModFreq )
                obj.fwdModel.modFreq = newModFreq;
            else
                error('ModFreq should be numeric.')
            end
        end
        
        function obj = set.nPhotons( obj, newNPhotons )
            if isnumeric( newNPhotons )
                obj.fwdModel.modFreq = newNPhotons;
            else
                error('nPhotons should be numeric.')
            end
        end
        
        function obj = set.nRepetitions( obj, newNReps )
            if isnumeric( newNReps )
                obj.fwdModel.modFeq = newNReps;
            else
                error('nRepetitions should be numeric.')
            end
        end
        
        function out = get.gpuId( obj )
            out = obj.fwdModel.gpuId;
        end
        
        function out = get.image( obj )
            out = obj.fwdModel.image;
        end
        
        function out = get.probe( obj )
            out = obj.fwdModel.probe;
        end
        
        function out = get.modFreq( obj )
            out = obj.fwdModel.modFreq;
        end
        
        function out = get.nPhotons( obj )
            out = obj.fwdModel.nPhotons;
        end
        
        function out = get.nRepetitions( obj )
            out = obj.fwdModel.nRepetitions;
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
                    mask = obj.image.volume(:) == iLayer;
                    ext = nirs.getExtinctions( obj.probe.lambda(iLambda) );
                    
                    Js.HbO(lst,mask(:)) = Js.HbO(lst,mask(:)) * ext(1);
                    Js.HbR(lst,mask(:)) = Js.HbR(lst,mask(:)) * ext(2);
                    Js.water(lst,mask(:)) = Js.water(lst,mask(:)) * ext(3);
                    Js.lipd(lst,mask(:)) = Js.lipid(lst,mask(:)) * ext(4);
                    Js.cytC(lst,mask(:)) = Js.cytC(lst,mask(:)) * ext(5);
%                     Js.mus(lst,mask(:)) =
%                     -3*obj.optTissue(iLayer).kappa^2;
                    Js.mus(lst,mask(:)) = -1/3/obj.optTissue(iLayer).mus^2 * J.kappa(lst,mask(:));
                    
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
                if length(obj.optTissue(i).mus) == 1
                    mus = obj.optTissue(i).mus*ones(size(lambda));
                else
                    mus = obj.optTissue(i).mus;
                end
                ri = obj.optTissue(i).ri;
                optProp(i) = nirs.OpticalProperties( mua, mus, ri, lambda );
            end
            
            obj.fwdModel.optProp = optProp;
        end
    end
    
end


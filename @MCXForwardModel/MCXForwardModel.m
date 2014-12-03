classdef MCXForwardModel < nirs.ForwardModel
    %MCXFWDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties%( Constant )
        gpuId = 1;
    end
    
    properties
        image; 
        probe;
        optProp;
        
        directory = [nirs.defaultTmp() filesep 'tmp'...
            filesep num2str(randi(2^32-1))];
        cleanup = true;

        modFreq = 110e6;
        nPhotons = 1e7;
        nRepetitions = 1;
        nTimeGates = 32;
        timeStep = 1/110e6/32;
    end
    
    properties(SetAccess = private)
        nLayers;
    end
    
    methods
        %% Constructor
        function obj = MCXForwardModel( varargin )
            if nargin > 0
                obj.image = varargin{1};
            end
            
            if nargin > 1
                obj.optProp = varargin{2};
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
                obj.image = newImage;
                obj.nLayers = max( obj.image.volume(:) );
            else
                error('Image should be of Image class with uint8 volume.')
            end
        end
        
        function obj = set.probe( obj, newProbe )
            if isa( newProbe,'nirs.Probe' )
                obj.probe = newProbe;
            else
                error('Probe should be of Probe class.')
            end
        end
        
        function obj = set.optProp( obj, newProp )
            if isa( newProp,'nirs.OpticalProperties' )
                obj.optProp = newProp;
            else
                error('Optical properties should be of OpticalProperties class.')
            end
        end
        
        function obj = set.directory( obj, newDir )
            if isstr( newDir )
                obj.directory = newDir;
                obj.cleanup = false;
            else
                error('Directory should be a string.')
            end
        end
        
        %% Methods
        meas = measurement( obj );
        [J,meas] = jacobian( obj );
        saveFluence( obj );
        
        function data = getMeasFromFluence(obj,pos,fluence)
            for iPos = 1:size(pos,1)
%                 x = pos(iPos,1) / obj.image.dim(1);
%                 y = pos(iPos,2) / obj.image.dim(1);
%                 z = pos(iPos,3) / obj.image.dim(1);
% 
%                 [X, Y, Z] = meshgrid( ...
%                     max(floor(x-1),1):min(ceil(x+1),size(fluence,1)),...
%                     max(floor(y-1),1):min(ceil(y+1),size(fluence,2)),...
%                     max(floor(z-1),1):min(ceil(z+1),size(fluence,3)) );
% 
%                 V = fluence( sub2ind( size(fluence), X, Y, Z ) );
%                 
%                 mask = obj.image.volume( sub2ind( size(obj.image.volume), X, Y, Z ) ) > 0;
%                 
%                 V = log(V(mask));
%                 V = exp( median(real(V)) + median(imag(V))*1i );
%                 
%                 data(iPos) = V;
                
%                 if length(V) == 1
%                     data(iPos) = V;
%                 else
%                     data(iPos) = interp3( X, Y, Z, V, x, y, z );
%                 end

x = pos(iPos,1) / obj.image.dim(1);
y = pos(iPos,2) / obj.image.dim(1);
z = pos(iPos,3) / obj.image.dim(1);

k = 3;

[X, Y, Z] = meshgrid( ...
    max(floor(x-k),1):min(ceil(x+k),size(fluence,1)),...
    max(floor(y-k),1):min(ceil(y+k),size(fluence,2)),...
    max(floor(z-k),1):min(ceil(z+k),size(fluence,3)) );

V = fluence( sub2ind( size(fluence), X, Y, Z ) );

Vr = real(log(V));
Vi = imag(log(V));

mask = obj.image.volume( sub2ind( size(obj.image.volume), X, Y, Z ) ) > 0;

ar = pinv([X(mask) Y(mask) Z(mask) ...
    X(mask).*Y(mask) X(mask).*Z(mask) Y(mask).*Z(mask) ...
    X(mask).^2 Y(mask).^2 Z(mask).^2 ones(size(Z(mask)))]) ...
    * Vr(mask);

ai = pinv([X(mask) Y(mask) Z(mask) ...
    X(mask).*Y(mask) X(mask).*Z(mask) Y(mask).*Z(mask) ...
    X(mask).^2 Y(mask).^2 Z(mask).^2 ones(size(Z(mask)))]) ...
    * Vi(mask);

thisR = [x y z x*y x*z y*z x^2 y^2 z^2 1]*ar;
thisI = [x y z x*y x*z y*z x^2 y^2 z^2 1]*ai;

data(iPos) = exp( thisR + thisI * 1i );
            end
        end
        
    end
    
    %% Private Methods
    methods( Access = private )
        cfg = getConfig( obj, idx, idxType );
%         out = getFluence( obj, srcField, detPos );
%         out = getDataFromPhotons( obj, detps,iLambda );
    end
    
end


classdef SlabForwardModel < nirs.ForwardModel
    %SLABFORWARDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        probe;
        optProp;
        modFreq = 110e6;
    end
    
    properties( Constant )
        nLayers = 1;
    end
    
    methods
        %% Constructor
        function obj = SlabForwardModel( varargin )
           if nargin > 0
               obj.probe = varargin{1};
           end
           
           if nargin > 1
               obj.optProp = varargin{2};
           end
           
           if nargin > 2
               obj.modFreq = varargin{3};
           end
           
           if nargin > 3
               error('Too many input arguments.')
           end
        end
        
        %% Set/Get
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
        
        %% Methods
        function meas = measurement( obj )
            d = obj.probe.distances;
            w = obj.modFreq*2*pi;

            for iLink = 1:size( obj.probe.link,1 )
                iLambda = obj.probe.link(iLink,3);
                k = obj.optProp.kappa( iLambda );
                v = obj.optProp.c( iLambda );
                mua = obj.optProp.mua( iLambda );
                
%                 logAC(1,iLink) = -log(d(iLink)^2) - d(iLink) * sqrt(mua/k) * (1 + (w/v/mua)^2)^.25 * cos( .5*atan(w/v/mua) );
%                 phi(1,iLink) = - d(iLink) * sqrt(mua/k) * (1 + (w/v/mua)^2)^.25 * sin( .5*atan(w/v/mua) );

                Vplus = sqrt(sqrt(1+(w/mua/v)^2) + 1);
                Vminus = sqrt(sqrt(1+(w/mua/v)^2) - 1);
                
                logAC(1,iLink) = - d(iLink) * sqrt(mua/2/k) * Vplus - log(d(iLink)^3) - ...
                    log(sqrt(1+d(iLink)*sqrt(2*mua/k)*Vplus + d(iLink)^2 * mua/k * sqrt(1 + (w/v/mua)^2)));
                phi(1,iLink) = -d(iLink) * sqrt(mua/2/k) * Vminus - atan(d(iLink)*sqrt(mua/2/k)*Vminus/(1+d(iLink)*sqrt(mua/2/k)*Vplus));
            end
            
            meas = nirs.Data( 0,exp(logAC + 1i*phi),obj.probe,obj.modFreq,'Analytical slab forward model.');
        end
        
        function [J,meas] = jacobian( obj )
            for iLambda = 1:length(obj.optProp.lambda)
                meas = obj.measurement();
                y0 = meas.data.';

%                 % ANALYTIC JACOBIAN IS LESS ACCURATE
%                 d = obj.probe.distances;
%                 w = obj.modFreq*2*pi;
%                 k = obj.optProp.kappa( iLambda );
%                 v = obj.optProp.c( iLambda );
%                 mua = obj.optProp.mua( iLambda );
%                 J(iLambda).mua = - d * ( w*sin(0.5*atan(w/mua/v)) + mua*v*cos(0.5*atan(w/mua/v)) ) ...
%                     / ( 2*k*mua*v*sqrt(mua/k)*(1 + (w/mua/v)^2)^0.75 )...
%                     - 1i * d * ( mua*v*sin(0.5*atan(w/mua/v)) - w*cos(0.5*atan(w/mua/v)) )...
%                     / ( 2*k*mua*v*sqrt(mua/k)*(1 + (w/mua/v)^2)^0.75 );
%                 J(iLambda).kappa = d * ( sqrt(mua/k)*(1 + (w/mua/v)^.25)*cos(0.5*atan(w/mua/v)) ) / 2/k ...
%                     + 1i * d * ( sqrt(mua/k)*(1 + (w/mua/v)^.25)*sin(0.5*atan(w/mua/v)) ) / 2/k;

                %dy/dmua
                thisModel = obj;
                thisModel.optProp.mua(iLambda) = 1.0001*thisModel.optProp.mua(iLambda);
                newMeas = thisModel.measurement();
                y1 = newMeas.data.';

                J(iLambda).mua = (y1 - y0)/(.0001*obj.optProp.mua(iLambda)) ./ y0;

                %dy/dk
                thisModel = obj;
                thisModel.optProp.kappa(iLambda) = 1.0001*thisModel.optProp.kappa(iLambda);
                newMeas = thisModel.measurement();
                y1 = newMeas.data.';

                J(iLambda).kappa = (y1 - y0)/(.0001*obj.optProp.kappa(iLambda)) ./ y0;
            end
        end
    end
    
    
    
end


classdef SlabForwardModel 
    %SLABFORWARDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        probe;
        prop;
        Fm = 110;
    end
    
    methods
        %% Constructor
        function obj = SlabForwardModel( probe, prop, Fm )
           if nargin > 0, obj.probe = probe; end
           if nargin > 1, obj.prop = prop; end
           if nargin > 2, obj.Fm = Fm; end
        end
        
        %% Methods
        function meas = measurement( obj )
            d = obj.probe.distances;
            w = obj.Fm*2*pi * 1e6;

            for iLink = 1:size( obj.probe.link,1 )
                iLambda = obj.probe.link(iLink,3);
                k = obj.prop.kappa( iLambda );
                v = obj.prop.v( iLambda );
                mua = obj.prop.mua( iLambda );

                Vp = sqrt(sqrt(1+(w/mua/v)^2) + 1);
                Vm = sqrt(sqrt(1+(w/mua/v)^2) - 1);
                
             	G = 1 + d(iLink) * sqrt(2*mua/k) * Vp + d(iLink).^2 * mua/k * sqrt(1+(w/v/mua)^2);
                
                logAC(1,iLink) = - d(iLink) * sqrt(mua/2/k) * Vp - log(d(iLink)^3) + log(sqrt(G));
              
                phi(1,iLink) = -d(iLink) * sqrt(mua/2/k) * Vm + atan(d(iLink)*sqrt(mua/2/k)*Vm/(1+d(iLink)*sqrt(mua/2/k)*Vp));
            end
            
            meas = nirs.Data( exp(logAC + 1i*phi),...
                obj.probe,...
                0, ...
                obj.Fm, ...
                'Analytical slab forward model.');
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


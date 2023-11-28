classdef KalmanAR < handle
    properties%( SetAccess = private )
        kf = nirs.realtime.util.KalmanFilter.empty; 	% kalman filters
        P;                      	% model orders
        Pmax;
        BIC;                        % information criterion
        resid;
    end
    
    properties( Dependent )
        a;
    end
    
    properties( Access = private )
        ybuffer;
        n = 0;
    end
        
    methods
        % boa constructer
        function obj = KalmanAR( P )
            % model orders
            obj.P = P;
            obj.Pmax = max(P);
            
            % create kalman filters
            for i = 1:length(P)
                obj.kf(i) = nirs.realtime.util.KalmanFilter( 0*eye(P(i)+1) );
                obj.kf(i).P = 10*eye(P(i)+1);
                obj.kf(i).B(2) = 0.5;
               
            end
            
            % allocate BIC array
            obj.BIC = zeros(length(P),1);
            
            obj.ybuffer = zeros(max(P)+1,1);
        end
 
        % update method
        function update(obj,y) % takes the last Pmax points

            obj.n = obj.n + 1;
            
            for i = 1:length(obj.kf)
                obj.kf(i).update(y,[1; obj.ybuffer(1:obj.P(i))]');
                obj.BIC(i) = obj.n * log( obj.kf(i).E.^2 ) + (obj.P(i)+1)*log( obj.n );
                obj.resid(i)=obj.kf(i).resid;
            end
            % 
            % for i = 1:length(obj.kf)
            %     obj.kf(i).update(y(1),[1; y(2:obj.P(i)+1)]');
            %     obj.BIC(i) = obj.n * log( obj.kf(i).E.^2 ) + (obj.P(i)+1)*log( obj.n );
            % 
            % end
            
            obj.ybuffer = [y; obj.ybuffer(1:end-1)];
 
        end

        function set.a(obj,a)
            if obj.n > 10*obj.Pmax
                [~,i] = min(obj.BIC);
                obj.kf(i).B=a;
            else
                obj.kf(end).B=a;
            end
        end
        
        function a = get.a(obj)
            if obj.n > 10*obj.Pmax
                [~,i] = min(obj.BIC);
                a = obj.kf(i).B;
            else
                a = obj.kf(end).B;
            end
        end
    end
    
end


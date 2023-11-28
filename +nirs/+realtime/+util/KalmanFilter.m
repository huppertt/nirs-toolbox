classdef KalmanFilter < handle
    properties
        Q = 0; 	% process noise
        P;      % covariance matrix
    end
    
    properties%( SetAccess = protected ) 
        resid; % current value of residual
        B;      % parameters
        E  = 0; % measurement noise
        Em = 0;
        n  = 0;	% number of updates
        np;     % number of parameters
    end
    
    properties( Access = private )
        S;
        K;
    end
    
    methods
        function obj = KalmanFilter( Q )
            if size(Q,1) ~= size(Q,2)
                error;
            end
            obj.Q = Q;
            obj.np = length(Q);
            
            obj.P = 100*eye(obj.np);
            obj.B = zeros(obj.np,1);
        end
        
        function update(obj,y,X,F,u)
            if size(y,1) ~= size(X,1)
                error;
            end
            
            if nargin < 4
                F = 1;
            end
            if nargin < 5
                u = 0;
            end
            
            obj.n = obj.n + 1;
            
        	% update measurement noise
            ehat = y - X*obj.B;
            obj.updateE( ehat );
            obj.resid=ehat;
            if obj.n > 0; %10
                % prediction step
                obj.prediction_step(F,u);

                % update B and P
                obj.update_step(ehat,X);
            end
        end
    end  
    methods( Access = protected )
        function prediction_step(obj,F,u)
            obj.B = F*obj.B + u*obj.B;
            obj.P = F*obj.P*F' + obj.Q;
        end

        function updateE( obj,ehat )
            %obj.Em = obj.Em*(obj.n-1)/obj.n + ehat/obj.n;
            obj.E  = obj.E *(obj.n-1)/obj.n + (1.253 * abs(ehat-obj.Em))/obj.n;
        end

        function update_step(obj,ehat,X)
            obj.S = X*obj.P*X' + obj.E*obj.E';
            %obj.K = obj.P*X'/obj.S;
            obj.K = obj.P*X'*pinv(obj.S);

            obj.B = obj.B + obj.K*ehat;
            obj.P = obj.P - obj.K*X*obj.P;  
        end
    end    
end
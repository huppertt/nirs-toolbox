classdef RecursiveLeastSquares < handle
    properties
        P;          % covariance matrix
        B;          % parameters
        lambda = 1; %
        
        E = 0;      % residual std dev
        Em = 0;     % mean of residual
        n = 0;      % count of updates
    end
    
    methods
        % constuctor
        function obj = RecursiveLeastSquares( P )
            if size(P,1) ~= size(P,2)
                error;
            end
            
            obj.P = P;
            obj.B = zeros(size(obj.P,1),1);
        end
        
        % update function 
        function update(obj,y,X)
            obj.n = obj.n + 1;
            
            % RLS steps
            K = obj.P*X' / (X*obj.P*X' + obj.lambda);
            obj.P = obj.P - K*X*obj.P;
            obj.B = obj.B + K*(y-X*obj.B);
            
            % statistics
            ehat = y - X*obj.B;
            obj.Em = obj.Em*(obj.n-1)/obj.n + ehat/obj.n;
            obj.E  = obj.E *(obj.n-1)/obj.n + (1.253 * abs(ehat-obj.Em))/obj.n;
        end
    end  
end
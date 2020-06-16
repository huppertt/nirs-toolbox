classdef MotionCorrect < handle
   properties
       Tune=1;
   end
   properties(Hidden=true)
        model=4;
        Q1=[];
        R1=[];
        P1=[];
        X1=[];
        I1=[];
        LagMtx=[];
        
        Q2=[];
        R2=[];
        P2=[];
        X2=[];
        I2=[];
   end
    
    methods
        function obj = MotionCorrect()
        end
        
        function d = update(obj,d,t)
            
            if(isempty(obj.Q1))
                obj.I1=speye(obj.model*size(d,2),obj.model*size(d,2));
                obj.X1=zeros(obj.model*size(d,2),1);
                obj.LagMtx=zeros(obj.model,size(d,2));
                obj.Q1 = obj.I1*0;
                obj.R1 = speye(size(d,2));
                obj.P1 = obj.I1;
                
                obj.I2=speye(size(d,2),size(d,2));
                obj.X2=d';
                obj.Q2 = obj.I2*.5;
                obj.R2 = 10*speye(size(d,2));
                obj.P2 = obj.I2*5;
                
            end
            
            % first update the AR model
            %for i=1:size(d,2)
            H=zeros(size(d,2),size(d,2)*obj.model);
            for i=1:size(obj.LagMtx,2)
                H(i,i:size(d,2):end)=obj.LagMtx(:,i);
            end
            yhat = H*obj.X1;
            innov = d'-yhat;
            obj.P1=obj.P1+obj.Q1;
            K = obj.P1*H'*inv(H*obj.P1*H'+obj.R1);
            obj.X1 = obj.X1 + K*innov;
            obj.P1 = (obj.I1-K*H)*obj.P1;
            
            obj.LagMtx(2:end,:)=obj.LagMtx(1:end-1,:);
            
            % now update the weighted innovations model
            
            innov = innov-obj.X2;
            obj.P2=obj.P2+obj.Q2;
            K=obj.P2*inv(obj.P2+obj.R2);
            obj.X2= obj.X2+K*innov;
            obj.P2=(obj.I2-K)*obj.P2;
            
            d = (obj.X2+ H*obj.X1)';
            obj.LagMtx(1,:)=d';
           
            
        end
        
    end
    
    
end
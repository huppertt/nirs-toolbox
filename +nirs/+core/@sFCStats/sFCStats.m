classdef sFCStats
    %% sFC - Holds stats info for a connectivity model
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     probe        - Probe object describing measurement geometry
    %     demographics - Dictionary containing demographics info
    %     R            - Correlation values
    %     Z            - Fisher Z-transform
    %     p            - p-values f
    %     dfe          - degrees of freedom
    %     
    %  Methods:
    %     
    %     draw        - draws beta or tstat values on top of probe geometry
    %     table       - returns a table of all stats (minus full covariance)
  
    properties
        type            % connectivity model from +nirs/+sFC/
        description     % description of data (e.g. filename)      
        probe           % Probe object describing measurement geometry
        demographics    % Dictionary containing demographics info
        conditions;
        % Results storage
        R               % correlation value (depends on model)
        dfe          	% degrees of freedom
       
    end
    
    properties ( Dependent = true )       
        p           % p-value (depends on model)
        Z           % Fisher-Z transform
    end
    
    
    methods
      
         function Z = get.Z( obj )
             Z = abs(.5*log((1+obj.R)./(1-obj.R))).*sign(obj.R);
             Z(find(Z>6))=6;  % Fix the R=1 -> inf;  tanh(6) ~1 so cut there to keep the scale
             Z(find(Z<-6))=-6;
         end
        
         function p = get.p(obj)
             for idx=1:length(obj.conditions)
                 t(:,:,idx)=obj.R(:,:,idx)./sqrt((1-obj.R(:,:,idx).^2)/(max(obj.dfe(idx))-2));
                 p(:,:,idx)=2*tcdf(-abs(t),max(obj.dfe(idx)));
             end
         end
        
        function out = table( obj )
            %% table - returns a table of the regression stats
            link=obj.probe.link;
            
            if(~iscellstr(link.type))
                link.type=arrayfun(@(x){num2str(x)},link.type);
            end
            
            [i,j]=meshgrid(1:height(link),1:height(link));
           
            sourceFrom=link.source(i);
            detectorFrom=link.detector(i);
            typeFrom=link.type(i);
            
            sourceTo=link.source(j);
            detectorTo=link.detector(j);
            typeTo=link.type(j);
            out=table;
            for i=1:length(obj.conditions)
                cond=repmat({obj.conditions{i}},length(sourceFrom(:)),1);
                out = [out; table(cond,[sourceFrom(:)],[detectorFrom(:)],{typeFrom{:}}',...
                    [sourceTo(:)],[detectorTo(:)],{typeTo{:}}',...
                    obj.R(:),obj.Z(:),obj.p(:),...
                    'VariableNames',{'condition','SourceOrigin','DetectorOrigin','TypeOrigin',...
                    'SourceDest','DetectorDest','TypeDest','R','Z','pvalue'})];
            end
        end
        
        draw( obj, vtype, vrange, thresh );     

    end
    
    methods (Access = protected)
        newNames = transformNames( obj, T );
    end
  
end
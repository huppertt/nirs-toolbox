classdef ConnectivityStats
    %% ConnectivitySTATS - Holds stats info for a connectivity model
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     probe        - Probe object describing measurement geometry
    %     demographics - Dictionary containing demographics info
    %     variables    - ?????
    %     Grangers     - Grangers coefficients
    %     Pearsons      - Pearson coefficients
    %     p            - p-values for Grangers model
    %     p_pearson    - p-values for Pearsons model
    %     dfe          - degrees of freedom
    %     
    %  Methods:
    %     
    %     draw        - draws beta or tstat values on top of probe geometry
    %     table       - returns a table of all stats (minus full covariance)
  
    properties
        type            % Correlation or Grangers
        description     % description of data (e.g. filename)      
        probe           % Probe object describing measurement geometry
        demographics    % Dictionary containing demographics info

        dfe1          	% degrees of freedom #1
        dfe2          	% degrees of freedom #2
       
        Grangers_F      % F-value on the grangers 
        Pearsons
        
    end
    
    properties ( Dependent = true )
        Grangers    % Grangers
        Z   % Z-transfrom of Pearsons
        p           % p-value on Grangers or Pearsons
    end
    
    
    methods
      
        function obj = GtoF(obj,G)
        % Converts Grnager's to F-test
           %G = log(sqrt(F*dfe1/dfe2+1))
           % (exp(G)^2-1)*dfe2/dfe1 =F
           obj.Grangers_F = (exp(G).^2-1).*obj.dfe2./obj.dfe1;
           
%           G=log(mad(rr)/mad(ru));
%           F=((mad(rr)^2-mad(ru)^2)/df1)/(mad(ru)^2/df2);

                 
           
        end
        function obj = ZtoR(obj,Z)
            obj.Pearsons=tanh(Z);
        end
         function z = get.Z( obj )
             if(strcmp(obj.type,'Correlation'))
                z = .5*log((1+obj.Pearsons)./(1-obj.Pearsons));
             else
                 z=[];
             end
         end
        
         
        
        
        function p = get.p(obj)
         %   p=fcdf(1./obj.Grangers_F, obj.dfe2, obj.dfe1);
            if(strcmp(obj.type,'Grangers'))
                p=1-fcdf(obj.Grangers_F, obj.dfe2, obj.dfe1,'upper');
                p(find(isnan(p)))=1;
            else
                t=obj.Pearsons./sqrt((1-obj.Pearsons.^2)/(max(obj.dfe2(:))-2));
                p=2*tcdf(-abs(t),max(obj.dfe2(:)));
            end
        end
        
        function Grangers = get.Grangers( obj )
            try
                Grangers=log(sqrt(obj.Grangers_F.*obj.dfe1./obj.dfe2+1));
                Grangers(find(isnan(Grangers)))=0;
            catch
                Grangers=[];
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
            
            if(strcmp(obj.type,'Grangers'))
                out = table([sourceFrom(:)],[detectorFrom(:)],{typeFrom{:}}',...
                    [sourceTo(:)],[detectorTo(:)],{typeTo{:}}',...
                    obj.Grangers(:),obj.Grangers_F(:),obj.p(:),...
                    'VariableNames',{'SourceOrigin','DetectorOrigin','TypeOrigin',...
                    'SourceDest','DetectorDest','TypeDest','Grangers','F','pvalue'});
            else
                out = table([sourceFrom(:)],[detectorFrom(:)],{typeFrom{:}}',...
                    [sourceTo(:)],[detectorTo(:)],{typeTo{:}}',...
                    obj.Pearsons(:),obj.Z(:),obj.p(:),...
                    'VariableNames',{'SourceOrigin','DetectorOrigin','TypeOrigin',...
                    'SourceDest','DetectorDest','TypeDest','Pearsons','Z','pvalue'});
            end
            
        end
        
        stats = ttest( obj, c, b, names );
        
        draw( obj, vtype, vrange, thresh );     

    end
    
    methods (Access = protected)
        newNames = transformNames( obj, T );
    end
  
end
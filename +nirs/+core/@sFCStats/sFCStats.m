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
    %     
    %  Methods:
    %     
    %     draw        - draws the correlation values
    %     table       - returns a table of all stats (minus full covariance)
    %     graph       - returns a graph object from the connectivity model
    
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
        q
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
                 p(:,:,idx)=2*tcdf(-abs(t(:,:,idx)),max(obj.dfe(idx)));
             end
         end
        
         function q = get.q(obj)
            p=obj.p;
            q=ones(size(p));
            lst=find(p~=1);  
            q(lst)=nirs.math.BenjaminiHochberg(p(lst));
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
                    reshape(obj.R(:,:,i),size(cond)),reshape(obj.Z(:,:,i),size(cond)),reshape(obj.p(:,:,i),size(cond)),...
                    'VariableNames',{'condition','SourceOrigin','DetectorOrigin','TypeOrigin',...
                    'SourceDest','DetectorDest','TypeDest','R','Z','pvalue'})];
            end
        end
        
        function grph=graph(obj,vtype,thresh,flip)
            % converts to a graph type object
            
            if(nargin<2)
                vtype='Z';
            end
            if(nargin<3)
                thresh='p<0.05';
            end
            if(~isempty(strfind(thresh,'p')))
                stat='p';
            else
                stat='q';
            end
            thres=str2num(thresh(strfind(thresh,'<')+1:end));
            
            if(isempty(strfind(vtype,':')))
                type=obj(1).probe.link.type(1);
            else
                type=vtype(strfind(vtype,':')+1:end);
                vtype=vtype(1:strfind(vtype,':')-1);
            end
            
            if(nargin<4 || isempty(flip))
                flip=[1 1];
            end
            
            for i=1:length(obj)
                lst=find(ismember(obj(i).probe.link.type,type));
                Adj=obj(i).(vtype);
                Adj=(Adj.*(obj(i).(stat)<thres));
                Adj=Adj.*(~eye(size(Adj)));
                Adj=Adj(lst,lst);
                grph(i)=nirs.core.Graph(Adj);
                grph(i).description=obj(i).description;
                grph(i).demographics=obj(i).demographics;
                grph(i).probe=obj(i).probe;
                
                for k=1:length(lst)
                    j=lst(k);
                    Name=['Src-' num2str(obj(i).probe.link.source(j)) ...
                        ':Det-' num2str(obj(i).probe.link.detector(j)) ...
                        ' ' obj(i).probe.link.type{j}];
                    X=.5*(obj(i).probe.srcPos(obj(i).probe.link.source(j),1)+...
                        obj(i).probe.detPos(obj(i).probe.link.detector(j),1));
                    Y=.5*(obj(i).probe.srcPos(obj(i).probe.link.source(j),2)+...
                        obj(i).probe.detPos(obj(i).probe.link.detector(j),2));
                    Z=.5*(obj(i).probe.srcPos(obj(i).probe.link.source(j),3)+...
                        obj(i).probe.detPos(obj(i).probe.link.detector(j),3));
                    grph(i).nodeInfo.label{k}=Name;
                    grph(i).nodeInfo.X(k)=X;
                    grph(i).nodeInfo.Y(k)=Y;
                    grph(i).nodeInfo.Z(k)=Z;
                end
                
                
            end
           
            
            
            
        end
        S = ttest(obj, c, b, names);
        h=draw( obj, vtype, vrange, thresh,flip );     
        printAll( obj, vtype, vrange, thresh, folder, ext , flip);
    end
    
    
    methods (Access = protected)
        newNames = transformNames( obj, T );
    end
  
end
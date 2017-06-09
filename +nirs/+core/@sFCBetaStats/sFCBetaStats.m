classdef sFCBetaStats
    %% sFC - Holds stats info for a connectivity model
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     probe        - Probe object describing measurement geometry
    %     demographics - Dictionary containing demographics info
    %     b            - Correlation values
    %     covb            - Fisher Z-transform
    %     dfe          - degrees of freedom
    %
    %     tstat        - (dependent) t-stats of beta
    %     p            - (dependent) p-values of beta
    %     q            - (dependent) q-values of beta (false discovery rate)
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
        beta               % correlation value (depends on model)
        covb
        dfe          	% degrees of freedom
    end
    
    properties ( Dependent = true )
        tstat
        p           % p-value (depends on model)
        q
    end
    
    
    methods
              
        % t statistic calculation
        function tstat = get.tstat( obj )
            tstat = zeros(size(obj.beta));
            for idx=1:length(obj.conditions)
                tstat(:,:,idx) = obj.beta(:,:,idx) ./ sqrt(obj.covb(:,:,idx,idx));
            end
        end
        
        % p value calculation
        function p = get.p( obj )
            t = obj.tstat;
            p = ones(size(t));
            for idx=1:length(obj.conditions)
                p(:,:,idx) = 2*tcdf(-abs(t(:,:,idx)), obj.dfe(idx));
            end
        end
        
         function q = get.q(obj)
            p=obj.p;
            q=ones(size(p));
            for idx=1:length(obj.conditions)
                pcond = p(:,:,idx);
                qcond = ones(size(pcond));
                if nansum(nansum(pcond-pcond')) == 0
                    mask = triu(true(size(pcond)));
                else
                    mask = logical(eye(size(pcond,1)));
                end
                pcond(mask) = 1;
                lst=find(pcond~=1);  
                qcond(lst)=nirs.math.BenjaminiHochberg(pcond(lst));
                qcondT = qcond';
                qcond(mask) = qcondT(mask);
                q(:,:,idx) = qcond;
            end
         end
        
        % critical value
        function out = getCritT( obj, s )
            %% getCritT - returns critical t-stat
            % 
            % Args:
            %     s - string specifying statistical significance
            %         (e.g. 'p < 0.05' or 'q < 0.1')
        
            s = strtrim( strsplit( s, '<' ) );
            
            if s{1} == 'p'
                pcrit = str2num(s{2});
                out = - tinv( abs(pcrit)/2, obj.dfe );
                
            elseif s{1} == 'q'
                t = abs(obj.tstat(:))';
                q = obj.q(:);
                
                [~,idx] = unique(q);
                
                out = interp1(q(idx), t(idx), str2num(s{2}));
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
                    reshape(obj.beta(:,:,i),size(cond)),reshape(obj.tstat(:,:,i),size(cond)),reshape(obj.p(:,:,i),size(cond)),...
                    'VariableNames',{'condition','SourceOrigin','DetectorOrigin','TypeOrigin',...
                    'SourceDest','DetectorDest','TypeDest','beta','tstat','pvalue'})];
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
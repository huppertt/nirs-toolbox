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
        variables       % table describing regression coefficients
        description     % description of data (e.g. filename)      
        probe           % Probe object describing measurement geometry
        demographics    % Dictionary containing demographics info
        conditions;
        % Results storage
        R               % correlation value (depends on model)
        dfe          	% degrees of freedom
    end
        
   properties ( Hidden = true )     
        ZstdErr            %holds error covariance    
        pvalue_fixed;
    end
    
    properties ( Dependent = true )       
        p           % p-value (depends on model)
        Z           % Fisher-Z transform
        q
        t               %holds t-statistic value
        ishyperscan
        ishypersymm
        numunique
    end
    properties (Hidden = true, Dependent = true)
        uniqueinds
    end
    methods
         function ishyperscan = get.ishyperscan(obj)
             ishyperscan = any(strcmpi(obj.probe.link.Properties.VariableNames,'hyperscan'));
         end
        
         function Z = get.Z( obj )
             Z = abs(.5*log((1+obj.R)./(1-obj.R))).*sign(obj.R);
             Z(Z>6)=6;  % Fix the R=1 -> inf;  tanh(6) ~1 so cut there to keep the scale
             Z(Z<-6)=-6;
         end
         
         function t = get.t(obj)

             dfe(1:length(obj.conditions)) = obj.dfe;
             
             for idx=1:length(obj.conditions)

                 n = dfe(idx);
                 
                 if(isempty(obj.ZstdErr))
                     t(:,:,idx)=obj.R(:,:,idx).*sqrt((n-2)./(1-obj.R(:,:,idx).^2));
                 else
                     t(:,:,idx)=obj.Z(:,:,idx)./sqrt(obj.ZstdErr(:,:,idx,idx));
                 end
             end
         end
         
         function p = get.p(obj)

             if(~isempty(obj.pvalue_fixed))
                if(isa(obj.pvalue_fixed,'nirs.bootstrapping.bootstrap_result'))
                    p=obj.pvalue_fixed.pvalue;
                    p=reshape(p,size(obj.t));
                else
                    p=obj.pvalue_fixed;  % hack to allow reuse in the Conn Seed model
                end
                return
             end
             
             for idx=1:length(obj.conditions)
                 
                 p(:,:,idx) = 2*nirs.math.tpvalue(-abs(obj.t(:,:,idx)),max(obj.dfe));
                 p(:,:,idx)=tril(p(:,:,idx),-1)+tril(p(:,:,idx),-1)'+eye(size(p(:,:,idx)));
                 
              
             end
         end
         
         function ishypersymm = get.ishypersymm(obj)
            if obj.ishyperscan
                hyper = obj.Z(1:end/2,end/2+1:end,:); % select the between-subject portion of matrix
                ishypersymm = (norm(reshape(hyper,size(hyper,1),[])-reshape(permute(hyper,[2 1 3]),size(hyper,1),[]))<eps(1)*10);
            else
                ishypersymm = false;
            end
         end
         
         function uniqueinds = get.uniqueinds(obj)
            p=obj.p;
            
            % Select upper triangle of whole matrix or of upper-right quadrant
            if(obj.ishypersymm)
                mask=repmat(triu(true(size(p,1)),size(p,1)/2),1,1,size(p,3));
            else
                mask=repmat(triu(true(size(p,1)),1),1,1,size(p,3));
            end
            
            % Mask out within-subject portion
            if(obj.ishyperscan)
                mask(1:end/2,1:end/2,:) = false;
                mask(end/2+1:end,end/2+1:end,:) = false;
            end
            
            uniqueinds=find(mask);
            
         end
         
         function numunique = get.numunique(obj)
             numunique = length(obj.uniqueinds);
         end
         
         function q = get.q(obj)
                        
            p=obj.p;
            q=ones(size(p));
            lst = obj.uniqueinds;
            
            q(lst)=nirs.math.BenjaminiHochberg(p(lst));
           
            if(obj.ishypersymm)
                for i=1:size(q,3)
                    hyper = q(1:end/2,end/2+1:end,i);
                    hyper = min(hyper,hyper');
                    q(1:end/2,end/2+1:end,i) = hyper;
                end
            end

            for i=1:size(q,3)
                q(:,:,i)=min(q(:,:,i),q(:,:,i)');
            end
            
         end
         
           function out = sorted( obj, colsToSortBy )
            %% sorted - returns sorted stats by columns in variables
            out = obj;
            if nargin < 2
                colsToSortBy = {'source', 'detector', 'type'};
            end
            
            if(length(obj)>1)
                for idx=1:length(obj)
                    out(idx)=sorted(obj(idx),colsToSortBy);
                end
                return
            end
            [out.probe.link, idx] = sortrows(out.probe.link, colsToSortBy);
            out.R = obj.R(idx,idx,:);
            if(~isempty(obj.ZstdErr))
                out.ZstdErr = obj.ZstdErr(idx, idx,:,:);
            end

            if(~isempty(out.pvalue_fixed))
                out.pvalue_fixed=out.pvalue_fixed(idx);
            end
            
            % Sort conditions
            [~,idx] = sort( out.conditions );
            out.conditions = out.conditions(idx);
            out.R = out.R(:,:,idx);
            if(~isempty(out.ZstdErr))
                out.ZstdErr = out.ZstdErr(:,:,idx,idx);
            end
            if length(out.dfe)==length(out.conditions)
                out.dfe = out.dfe(idx);
            end

            

        end
         
         
        function out = table( obj )
            %% table - returns a table of the regression stats
               
            if(length(obj)>1)
                out=[];
                cnt=[];
                for id=1:length(obj)
                    tt=obj(id).table;
                    out=nirs.util.safe_table_vcat(out,tt);
                    cnt(id)=height(tt);
                end
                out2=[];
                for id=1:length(obj)
                    out2=nirs.util.safe_table_vcat(out2,repmat(nirs.createDemographicsTable(obj(id)),cnt(id),1));
                end
                out=[out2 out];
                return;
            end

            link=obj.probe.link;
            
            if(~iscellstr(link.type))
                link.type=arrayfun(@(x){num2str(x)},link.type);
            end
            
            [i,j]=meshgrid(1:height(link),1:height(link));
           
            if(isa(obj.probe,'nirs.core.ProbeROI') || ismember('ROI',link.Properties.VariableNames))
                ROIFrom=link.ROI(i);
                typeFrom=link.type(i);
                
                ROITo=link.ROI(j);
                typeTo=link.type(j);
                out=table;
                for i=1:length(obj.conditions)
                    cond=repmat({obj.conditions{i}},length(ROIFrom(:)),1);
                    out = [out; table(cond,[ROIFrom(:)],{typeFrom{:}}',...
                        [ROITo(:)],{typeTo{:}}',...
                        reshape(obj.R(:,:,i),size(cond)),reshape(obj.Z(:,:,i),size(cond)),...
                        reshape(obj.t(:,:,i),size(cond)),reshape(obj.p(:,:,i),size(cond)),reshape(obj.q(:,:,i),size(cond)),...
                        'VariableNames',{'condition','ROIOrigin','TypeOrigin',...
                        'ROIDest','TypeDest','R','Z','t','pvalue','qvalue'})];
                end
            else
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
                        reshape(obj.R(:,:,i),size(cond)),reshape(obj.Z(:,:,i),size(cond)),...
                        reshape(obj.t(:,:,i),size(cond)),reshape(obj.p(:,:,i),size(cond)),reshape(obj.q(:,:,i),size(cond)),...
                        'VariableNames',{'condition','SourceOrigin','DetectorOrigin','TypeOrigin',...
                        'SourceDest','DetectorDest','TypeDest','R','Z','t','pvalue','qvalue'})];
                end
            end
        end
        
        function seedGUI(obj)
            % wrapper to call the seed based connectivity GUI
        end
        
        function seed(obj,channel, vtype, vrange, thresh)
            % draws the stats map of regions connected to a seed channel
            
            % allows R or Z
            if(nargin<3 || isempty(vtype))
                vtype='R';
            end
            
            if(nargin < 5)
                thresh='p<0.05';
            end
            if(nargin<4)
                vrange=[];
            end
          
            
            if(strcmp(lower(vtype),'r'))
                values=obj.R(channel,:)';
            elseif(strcmp(lower(vtype),'z'))
                 values=obj.Z(channel,:)';
            else
                warning('unknown type.  Must be R or Z');
                return;
            end
                
            
            cond=repmat({['Channel_' num2str(channel)]},height(obj.probe.link),1);
            
            c=nirs.core.ChannelStats;
            c.probe=obj.probe;
            c.beta=values;
            c.variables=[c.probe.link table(values,cond,'VariableNames',{'beta','cond'})];
            c.pvalue_fixed=obj.p(channel,:)';
            c.covb=eye(length(values));
            hh=c.draw('beta', vrange, thresh);
            
            types=unique(obj.probe.link.type);
            
            
            for i=1:length(hh)
                
                
                if(obj.probe.link(channel,:).type==types(i))
                    figure(hh(i));
                hold on;
                sI=obj.probe.link.source(channel);
                 dI=obj.probe.link.detector(channel);
                src=obj.probe.srcPos(sI,:);
                det=obj.probe.detPos(dI,:);
                ss(1)=scatter(src(1),src(2));
                ss(2)=scatter(det(1),det(2));
                 
                set(ss,'SizeData',580);
                set(ss,'MarkerFaceColor','k');
                    h=line([src(1) det(1)],[src(2) det(2)],'Linewidth',10,'Color','k');
                end
                
                hold off;
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
                Adj(logical(eye(size(Adj))))=1;
                grph(i)=nirs.core.Graph(Adj);
                grph(i).description=obj(i).description;
                grph(i).demographics=obj(i).demographics;
                grph(i).probe=obj(i).probe;
                
                for k=1:length(lst)
                    j=lst(k);
                    if(isa(obj(i).probe,'nirs.core.ProbeROI'))
                         Name=['ROI-' obj(i).probe.link.ROI{j} ...
                        ' ' obj(i).probe.link.type{j}];
                        X=k*10;
                        Y=0;
                        Z=0;
                    
                    else
                    Name=['Src-' num2str(obj(i).probe.link.source(j)) ...
                        ':Det-' num2str(obj(i).probe.link.detector(j)) ...
                        ' ' obj(i).probe.link.type{j}];
                    X=.5*(obj(i).probe.srcPos(obj(i).probe.link.source(j),1)+...
                        obj(i).probe.detPos(obj(i).probe.link.detector(j),1));
                    Y=.5*(obj(i).probe.srcPos(obj(i).probe.link.source(j),2)+...
                        obj(i).probe.detPos(obj(i).probe.link.detector(j),2));
                    Z=.5*(obj(i).probe.srcPos(obj(i).probe.link.source(j),3)+...
                        obj(i).probe.detPos(obj(i).probe.link.detector(j),3));
                    
                    end
                   
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
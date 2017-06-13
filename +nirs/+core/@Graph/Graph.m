classdef Graph
    % This class holds node/edge values information
    %
    % Properties:
    %     description  - description of data (e.g. filename)
    %     demographics - Dictionary containing demographics info
    %     nodeInfo     - table with node positions and labels
    %     adjacency    -
    %     edgeInfo     -table with edge info
    %
    %  Methods:
    %     draw
    %     pagerank
    
    properties
        description;  % description of data (e.g. filename)
        demographics; % Dictionary containing demographics info
        nodeInfo;     % table with node positions and labels
        adjacency;    %
        edgeInfo;     %table with edge info
        probe;
    end
    methods
        function obj = Graph(Adj,n,label)
            % Constructor function to create Graph object
            % inputs:
            %   n - position of nodes
            %   Adj - adjacency matrix
            %   labels - labels of the edges
            
            if(nargin>0)
                obj.adjacency=Adj;
            else
                return
            end
            if(nargin<2 || isempty(label))
                theta=linspace(0,2*pi,size(Adj,1)+1)';
                [n(:,1),n(:,2)]=pol2cart(theta(1:end-1),ones(size(Adj,1),1));
                n(:,3)=0;
            end
            if(nargin<3 || isempty(label))
                for i=1:size(n,1)
                    label{i,1}=['node-' num2str(i)];
                end
            end
            
            X=n(:,1);
            Y=n(:,2);
            Z=n(:,3);
            value=ones(size(n,1),1);
            obj.nodeInfo=table(label,X,Y,Z,value);
            
            obj.edgeInfo=table;
            
            obj=obj.weightededge;
        end
        
        function obj = weightededge(obj)
            cnt=1;
            orig=[]; value=[]; dest=[]; directional=[];
            obj.nodeInfo.value(:)=1;
            for i=1:size(obj.adjacency,1)
                for j=1:size(obj.adjacency,2)
                    if(i~=j)
                        if(obj.adjacency(i,j)~=0)
                            orig(cnt,1)=i;
                            dest(cnt,1)=j;
                            value(cnt,1)=obj.adjacency(i,j);
                            if(obj.adjacency(i,j)==obj.adjacency(j,i))
                                directional(cnt,1)=false;
                            else
                                directional(cnt,1)=true;
                            end
                            cnt=cnt+1;
                            
                        end
                    end
                end
            end
            obj.edgeInfo=table(value,orig,dest,directional);
        end
        function obj = pagerank(obj,d,falff)
            if(nargin<2 || isempty(d))
                d=.85;
            end
            if(nargin<3 || isempty(falff))
                N = size(obj.adjacency,1);
                falff = ones(N,1)/N;
            end
            pgrank=pagerank_centrality(obj.adjacency,d,falff);
            obj.nodeInfo.value=pgrank;
            obj.edgeInfo.value(:)=1;
            obj.edgeInfo.directional(:)=false;
            
        end
        function obj = edge_betweenness(obj)
            [e,b]=edge_betweenness_bin(obj.adjacency);
            obj.nodeInfo.value=b;
            obj.edgeInfo.value=e(find(obj.adjacency~=0));
        end
        
        function obj = degrees(obj)
            a=degrees_dir(obj.adjacency);
            obj.nodeInfo.value=a(:);
            obj.edgeInfo.value(:)=1;
            obj.edgeInfo.directional(:)=false;
        end
        
        
        function obj=efficiency(obj)
            a=efficiency_wei(obj.adjacency,1);
            obj.nodeInfo.value=a;
            obj.edgeInfo.value(:)=1;
            obj.edgeInfo.directional(:)=false;
        end
        
        function h=draw(obj)
            
            if(~isempty(obj.probe))
                l=obj.probe.draw;
                if(isa(l(1),'matlab.graphics.chart.primitive.Scatter'))
                    set(l,'MarkerFaceColor',[.85 .85 .85],'MarkerEdgeColor',[.85 .85 .85])
                else
                    set(l,'color',[.85 .85 .85])
                end
            end
            
            hold on;
            caxis([-1 1]);
            c=max(abs(obj.edgeInfo.value));
            cmap=colormap(jet);
            
            if(length(unique(obj.edgeInfo.value))==1)
                cmap(:)=0;
            end
            
            cIdx=linspace(-c,c,size(cmap,1));
            
            for i=1:height(obj.edgeInfo)
                for j=1:3; col(j) = interp1(cIdx,cmap(:,j),obj.edgeInfo.value(i)); end;
                pos1=[obj.nodeInfo.X(obj.edgeInfo.orig(i))...
                    obj.nodeInfo.Y(obj.edgeInfo.orig(i))...
                    obj.nodeInfo.Z(obj.edgeInfo.orig(i))];
                pos2=[obj.nodeInfo.X(obj.edgeInfo.dest(i))...
                    obj.nodeInfo.Y(obj.edgeInfo.dest(i))...
                    obj.nodeInfo.Z(obj.edgeInfo.dest(i))];
                
                if(obj.edgeInfo.directional(i))
                    % todo: draw an arrow
                else
                    l=line([pos1(1) pos2(1)],[pos1(2) pos2(2)],[pos1(3) pos2(3)]);
                    set(l,'color',col,'linewidth',1);
                end
            end
            if(length(unique(obj.edgeInfo.value))>1)
                cb1=colorbar('WestOutside');
                
                if(all(obj.edgeInfo.value>=0))
                    set(cb1,'Limits',[0 1]);
                    set(cb1,'TickLabels',num2str(round(linspace(0,c,length(get(cb1,'TickLabels')))'*100)/100));
                else
                    set(cb1,'TickLabels',num2str(round(linspace(-c,c,length(get(cb1,'TickLabels')))'*100)/100));
                    
                end
            end
            
            c=max(abs(obj.nodeInfo.value));
            c=.3;
            cmap=colormap(jet);
            if(length(unique(obj.nodeInfo.value))==1)
                cmap(:)=0;
            end
            
            cIdx=linspace(-c,c,size(cmap,1));
            hold on;
            for i=1:height(obj.nodeInfo)
                for j=1:3; col(j) = interp1(cIdx,cmap(:,j),obj.nodeInfo.value(i)); end;
                
                s=scatter3(obj.nodeInfo.X(i),obj.nodeInfo.Y(i),obj.nodeInfo.Z(i),64,col,'filled');
                
            end
            
            if(length(unique(obj.nodeInfo.value))>1)
                cb2=colorbar('EastOutside');
                
                if(all(obj.nodeInfo.value>=0))
                    set(cb2,'Limits',[0 1]);
                    set(cb2,'TickLabels',num2str(round(linspace(0,c,length(get(cb2,'TickLabels')))'*100)/100));
                else
                    set(cb2,'TickLabels',num2str(round(linspace(-c,c,length(get(cb2,'TickLabels')))'*100)/100));
                
                end
                
            end
            axis off;
            
            hold off;
            
        end
        
    end
end

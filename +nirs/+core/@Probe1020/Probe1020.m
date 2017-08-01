classdef Probe1020 < nirs.core.Probe
    %% Probe1020.  This is a class derived from the core Probe class which adds
    % methods for registration, display, and co-registration based on the 10-20
    % landmarking
    
    properties
        optodes_registered;   % Registered version of optodes
        braindepth;  % depth of brain for modeling
        opticalproperties;  % optical properties for modeling
        
        defaultdrawfcn; % Drawing function to use
        
    end
    properties( Dependent = true )
        headcircum;  % circumference of head  % make dependent
        AP_arclength; % Anterior-Posterior arc length
        LR_arclength; % Left-right arclength
        
    end
    properties (Access = private)
        zoom; % Flag to zoom in or show full 10-20 probe
        
        AP_distance;  % Distance from Oz to Fpz
        LR_distance;  % Distance from LPA to RPA
        IS_distance;  % Distance from center to Cz
        labels;  % labels of 10-20 points
        pts1020;  % 10-20 points
        mesh;
        
    end
    
    methods
        
        function obj = Probe1020(probe,headsize)
            
            if(nargin>1)
                tbl=nirs.util.register_headsize(headsize);
            else
                tbl=nirs.util.list_1020pts('?');
            end
            obj.pts1020=[tbl.X tbl.Y tbl.Z];
            obj.labels=tbl.Name;
            
            % Find the default arc lengths
            pt(1,:)=obj.pts1020(find(ismember(lower(obj.labels),'lpa')),:);
            pt(2,:)=obj.pts1020(find(ismember(lower(obj.labels),'rpa')),:);
            pt(3,:)=obj.pts1020(find(ismember(lower(obj.labels),'cz')),:);
            pt(4,:)=obj.pts1020(find(ismember(lower(obj.labels),'nas')),:);
            pt(5,:)=obj.pts1020(find(ismember(lower(obj.labels),'iz')),:);
            
            obj.AP_distance=norm(pt(4,:)-pt(5,:));
            obj.LR_distance=norm(pt(1,:)-pt(2,:));
            obj.IS_distance=norm(pt(3,:)-.5*(pt(1,:)-pt(2,:)));
            
            obj.braindepth = 10;  % Fix and cite
            
            if(nargin>0 && ~isempty(probe))
                obj=obj.registerprobe(probe);
                lambda=unique(obj.optodes.type);
            else
                lambda=[808];
            end
            
            obj=obj.set_mesh;
            
            obj.opticalproperties = nirs.media.tissues.brain(0.7, 50,lambda);
            
            obj.defaultdrawfcn='draw1020';
        end
        
        function obj = regsister_mesh2probe(obj,mesh)
            % This reshapes the mesh to fit the current probe

            tbl=table(obj.labels,obj.pts1020(:,1),obj.pts1020(:,2),obj.pts1020(:,3),...
                'VariableNames',{'Name','X','Y','Z'});
            
            T = nirs.registration.cp2tform(mesh(1).fiducials,tbl);
            for idx=1:length(mesh)
                n=mesh(idx).nodes;
                n(:,4)=1;
                n=n*T;
                mesh(idx).nodes=n(:,1:3);
                
                if(~isempty(mesh(idx).fiducials))
                    p=[mesh(idx).fiducials.X mesh(idx).fiducials.Y mesh(idx).fiducials.Z];
                    p(:,4)=1;
                    p=p*T;
                    mesh(idx).fiducials.X=p(:,1);
                    mesh(idx).fiducials.Y=p(:,2);
                    mesh(idx).fiducials.Z=p(:,3);
                    
                    
                end
                
            end
            obj.mesh=mesh;
            
            
        end
        
              
        function headsize=get_headsize(obj)
            headsize=Dictionary();
            headsize('lpa-cz-rpa')=obj.LR_arclength;
            headsize('Iz-cz-nas')=obj.AP_arclength;
            headsize('circumference')=obj.headcircum;
            
        end
        function varargout=draw(obj,varargin)
            
            if(ismember('hyperscan',obj.link.Properties.VariableNames))
                p=obj;
                p.link=p.link(ismember(p.link.hyperscan,'A'),:);

                S={}; D={};
                for i=1:height(p.link)
                    if iscell(p.link.source(i))
                        source = p.link.source{i};
                        detector = p.link.detector{i};
                    else
                        source = p.link.source(i);
                        detector = p.link.detector(i);
                    end
                    for j=1:length(source)
                        s=['000' num2str(source(j))];
                        S{end+1}=['Source-' s(end-3:end)];
                        d=['000' num2str(detector(j))];
                        D{end+1}=['Detector-' d(end-3:end)];
                    end
                end

                p.optodes_registered=p.optodes_registered(ismember(p.optodes.Name,{S{:} D{:}}),:);
                p.optodes=p.optodes(ismember(p.optodes.Name,{S{:} D{:}}),:);
                obj = p;
            end
            
            if(~isempty(strfind(obj.defaultdrawfcn,'zoom')))
                obj.zoom=true;
            else
                obj.zoom=false;
            end
            
            if(~isempty(strfind(obj.defaultdrawfcn,'10-20 map')));
                l=draw1020interp(obj,varargin{:});
            elseif(~isempty(strfind(obj.defaultdrawfcn,'3D')));
               
                if(strfind(obj.defaultdrawfcn,'('))
                    v=obj.defaultdrawfcn(strfind(obj.defaultdrawfcn,'(')+1:end);
                    v=v(1:strfind(v,')')-1);
                else
                    v='frontal';
                end
                
                l=draw3d(obj,varargin{:});
                if(~isempty(strfind(obj.defaultdrawfcn,'mesh')))
                    mesh=obj.getmesh;
                    h=mesh.draw;
                    axis tight;
                    nirs.util.rotateview(get(h(1),'Parent'),v)
                 end
                
            elseif(~isempty(strfind(obj.defaultdrawfcn,'10-20')));
                l=draw1020(obj,varargin{:});
            else
                l=draw@nirs.core.Probe(obj,varargin{:});
            end
            if(nargout>0)
                varargout{1}=l;
            end
        end
        
        function str = get.defaultdrawfcn(obj)
            str = obj.defaultdrawfcn;
        end
        function obj = set.defaultdrawfcn(obj,str)
            
            if(nargin==1)
                return
            end
            
            allowed={'10-20','10-20 mercator projection map';...
                '10-20 zoom', '10-20 mercator with restricted view';...
                '10-20 map', '10-20 mercator with underlain image';...
                '10-20 map zoom', '10-20 mercator with underlain image';...
                '3D', '3D line drawing ';...
                '3D mesh', '3D line drawing overlain on mesh';...
                '3D mesh (frontal)', '3D line drawing overlain on mesh';...
                '3D mesh (left)', '3D line drawing overlain on mesh';...
                '3D mesh (right)', '3D line drawing overlain on mesh';...
                '3D mesh (superior)', '3D line drawing overlain on mesh';...
                '3D mesh (top)', '3D line drawing overlain on mesh';...
                '3D mesh (posterior)', '3D line drawing overlain on mesh';...
                '3D mesh (occipital)', '3D line drawing overlain on mesh';...
                '3D mesh (back)', '3D line drawing overlain on mesh';...
                '3D mesh (inferior)', '3D line drawing overlain on mesh';...
                '3D mesh (bottom)', '3D line drawing overlain on mesh';...
                '2D', '2D probe layout'};
            
            if(~isempty(str))
                idx=find(ismember(lower({allowed{:,1}}),lower(str)));
            else
                idx=[];
            end
            
            if(~isempty(idx))
                obj.defaultdrawfcn=allowed{idx,1};
            elseif(isempty(str) || strcmp(str,'?'))
                disp('Here are the options for drawing configurations');
                disp(allowed);
            end
        end
        
        
        function obj=swap_reg(obj)
            op=obj.optodes;
            obj.optodes=obj.optodes_registered;
            obj.optodes_registered=op;
        end
        function mesh = getmesh(obj);
            % create the spherical mesh
            mesh=obj.mesh;
        end
        
        function obj=set_mesh(obj);
            mesh(1) = nirs.util.spheremesh;
            mesh(2) = nirs.util.spheremesh;
            mesh(3) = nirs.util.spheremesh;
            
            pt(1,:)=obj.pts1020(find(ismember(obj.labels,'nas')),:);
            pt(2,:)=obj.pts1020(find(ismember(obj.labels,'Iz')),:);
            com=.5*sum(pt,1);
            
            fract=.2;
            mesh(1).nodes(:,1)=mesh(1).nodes(:,1)*obj.LR_distance/2+com(1);
            mesh(1).nodes(:,2)=mesh(1).nodes(:,2)*obj.AP_distance/2+com(2);
            mesh(1).nodes(:,3)=mesh(1).nodes(:,3)*obj.IS_distance+com(3);
            
            mesh(2).nodes(:,1)=mesh(2).nodes(:,1)*(obj.LR_distance/2-obj.braindepth*fract)+com(1);
            mesh(2).nodes(:,2)=mesh(2).nodes(:,2)*(obj.AP_distance/2-obj.braindepth*fract)+com(2);
            mesh(2).nodes(:,3)=mesh(2).nodes(:,3)*(obj.IS_distance-obj.braindepth*fract)+com(3);
            
            mesh(3).nodes(:,1)=mesh(3).nodes(:,1)*(obj.LR_distance/2-obj.braindepth)+com(1);
            mesh(3).nodes(:,2)=mesh(3).nodes(:,2)*(obj.AP_distance/2-obj.braindepth)+com(2);
            mesh(3).nodes(:,3)=mesh(3).nodes(:,3)*(obj.IS_distance-obj.braindepth)+com(3);
            
            
            [TR, TT] = icp(mesh(1).nodes',obj.pts1020');
            mesh(1).nodes=(TR'*mesh(1).nodes'-TT*ones(1,size(mesh(1).nodes,1)))';
            mesh(2).nodes=(TR'*mesh(2).nodes'-TT*ones(1,size(mesh(2).nodes,1)))';
            mesh(3).nodes=(TR'*mesh(3).nodes'-TT*ones(1,size(mesh(3).nodes,1)))';
            
            fidtbl=table(obj.labels,obj.pts1020(:,1),obj.pts1020(:,2),obj.pts1020(:,3),...
                repmat({'10-20'},length(obj.labels),1),...
                repmat({'mm'},length(obj.labels),1),...
                repmat(true,length(obj.labels),1),...
                'VariableNames',{'Name','X','Y','Z','Units','Type','Draw'});
            mesh(1).fiducials=fidtbl;
            mesh(1).transparency=.2;
            mesh(2).transparency=.1;
            obj.mesh=mesh;
            
        end
        
        function varargout=draw3d(obj,colors, lineStyles, axis_handle)
            link = obj.link(strcmp(obj.link.type,obj.link.type(1)),1:2);
            
            
            n = height(link);
            
            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end
            
            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            elseif size(lineStyles, 1) == 1
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            end
            
            if nargin < 4
                axis_handle = axes();
            end
            
            % Points from the probe
            Pos(:,1)=obj.optodes_registered.X;
            Pos(:,2)=obj.optodes_registered.Y;
            Pos(:,3)=obj.optodes_registered.Z;
            
            hold on;
            lstS=find(ismember(obj.optodes.Type,'Source'));
            scatter3(Pos(lstS,1),Pos(lstS,2),Pos(lstS,3),'filled','MarkerFaceColor','r')
            lstD=find(ismember(obj.optodes.Type,'Detector'));
            scatter3(Pos(lstD,1),Pos(lstD,2),Pos(lstD,3),'filled','MarkerFaceColor','b')
            
            for i=1:height(link)
                if iscell(link.source(i))
                    source = link.source{i};
                    detector = link.detector{i};
                else
                    source = link.source(i);
                    detector = link.detector(i);
                end
                for j=1:length(source)
                    s = source(j);
                    d = detector(j);
                    h(i)=line(Pos([lstS(s) lstD(d)],1),Pos([lstS(s) lstD(d)],2),Pos([lstS(s) lstD(d)],3),'Color', colors(i, :), lineStyles{i, :});
                    set(h(i),'UserData',[s d]);
                end
            end
            
            axis equal;
            
            if(nargout>0)
                varargout{1}=h;
            end
            
        end
        
        function varargout=draw1020interp(obj,colors, lineStyles, axis_handle)
            link = obj.link(strcmp(obj.link.type,obj.link.type(1)),1:2);
            n = height(link);
            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end
            
            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            elseif size(lineStyles, 1) == 1
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            end
            
            if nargin < 4
                axis_handle = axes();
            end
            
            h=draw1020(obj,[],[],axis_handle);
            delete(h);
            lstS=find(ismember(obj.optodes.Type,'Source'));
            lstD=find(ismember(obj.optodes.Type,'Detector'));
            
            lst=find(~ismember(obj.optodes_registered.Type,{'FID-anchor','FID-attractor'}));
            Pos(:,1)=obj.optodes_registered.X(lst);
            Pos(:,2)=obj.optodes_registered.Y(lst);
            Pos(:,3)=obj.optodes_registered.Z(lst);
            
            [x,y]=obj.convert2d(Pos);
            xylink = [];
            xycolors = [];
            for i=1:height(link)
                if iscell(link.source(i))
                    source = link.source{i};
                    detector = link.detector{i};
                else
                    source = link.source(i);
                    detector = link.detector(i);
                end
                for j = 1:length(source)
                    s=source(j);
                    d=detector(j);
                    xylink(end+1,1)=mean(x([lstS(s) lstD(d)]));
                    xylink(end,2)=mean(y([lstS(s) lstD(d)]));
                    xycolors(end+1,:) = colors(i,:);
                end
            end
            
            xlim=get(axis_handle,'XLim');
            ylim=get(axis_handle,'YLim');
            [x2,y2]=meshgrid(xlim(1):xlim(2),ylim(1):ylim(2));
            
            F = scatteredInterpolant(xylink(:,1),xylink(:,2),xycolors(:,1),'nearest','none');
            cm(:,:,1)=reshape(F(x2(:),y2(:)),size(x2));
            F.Values=xycolors(:,2);
            cm(:,:,2)=reshape(F(x2(:),y2(:)),size(x2));
            F.Values=xycolors(:,3);
            cm(:,:,3)=reshape(F(x2(:),y2(:)),size(x2));
            
            
            k=dsearchn(xycolors,reshape(cm,[],3));
            k=reshape(k,size(cm,1),size(cm,2));
            lst=find(ismember(k,find(strcmp(lineStyles(:,2),'--'))));
            
            mask=ones(size(cm,1),size(cm,2));
            mask(lst)=NaN;
            cm(:,:,1)=mask.*cm(:,:,1);
            cm(:,:,2)=mask.*cm(:,:,2);
            cm(:,:,3)=mask.*cm(:,:,3);
            
            i=imagesc(x2(:),y2(:),cm);
            h=draw1020(obj,[],[],axis_handle);
            set(h,'LineWidth',.1,'color',[.3 .3 .3]);
            set(i,'alphaData',~isnan(cm(:,:,1)));
            
            if(obj.zoom)
                dx=(max(x)-min(x))/20;
                dy=(max(y)-min(y))/20;
                set(gca,'Xlim',[min(x)-dx max(x)+dx]);
                set(gca,'Ylim',[min(y)-dy max(y)+dy]);
            end
            
            if(nargout>0)
                varargout{1}=h;
            end
        end
        
        function varargout=draw1020(obj,colors, lineStyles, axis_handle)
            % Code to draw the probe in 10-20 space
            
            if(~isempty(obj.link))
                link = obj.link(strcmp(obj.link.type,obj.link.type(1)),1:2);
            else
                link=[];
            end
            
            n = height(link);
            
            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end
            
            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 3}, [n 1]);
            elseif size(lineStyles, 1) == 1
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 3}, [n 1]);
            end
            
            if nargin < 4
                axis_handle = axes();
            end
            
            hold on;
            [x,y]=obj.convert2d(obj.pts1020);
            dx=-x(find(ismember(obj.labels,'Cz')));
            dy=-y(find(ismember(obj.labels,'Cz')));
            scatter(x+dx,y+dy,'filled','MarkerFaceColor',[.8 .8 .8]);
         
            
            
            
            % Todo-  draw the probe too
            
            if(~isempty(obj.optodes_registered))
                % Points from the probe
                lst=find(~ismember(obj.optodes_registered.Type,{'FID-anchor','FID-attractor'}));
                Pos(:,1)=obj.optodes_registered.X(lst);
                Pos(:,2)=obj.optodes_registered.Y(lst);
                Pos(:,3)=obj.optodes_registered.Z(lst);
                [x,y]=obj.convert2d(Pos);
                xop=x; yop=y;
                lstS=find(ismember(obj.optodes.Type,'Source'));
                scatter(x(lstS)+dx,y(lstS)+dy,'filled','MarkerFaceColor','r')
                lstD=find(ismember(obj.optodes.Type,'Detector'));
                scatter(x(lstD)+dx,y(lstD)+dy,'filled','MarkerFaceColor','b')
                
                for i=1:height(link)
                    if iscell(link.source(i))
                        source = link.source{i};
                        detector = link.detector{i};
                    else
                        source = link.source(i);
                        detector = link.detector(i);
                    end
                    for j=1:length(source)
                        s = source(j);
                        d = detector(j);
                        h(i)=line(x([lstS(s) lstD(d)])+dx,y([lstS(s) lstD(d)])+dy,'Color', colors(i, :), lineStyles{i, :});
                        set(h(i),'UserData',[s d]);
                    end
                end
                %%
            else
                h=[];
            end
            
            axis tight;
            axis equal;
            axis off;
            [x,y]=obj.convert2d(obj.pts1020);
            headradius=obj.headcircum/(2*pi);
            
           headradius=norm([x(find(ismember(obj.labels,'Fpz'))) y(find(ismember(obj.labels,'Fpz')))]-...
                [x(find(ismember(obj.labels,'Cz'))) y(find(ismember(obj.labels,'Cz')))]);
            
%              headradius=norm(obj.pts1020(find(ismember(obj.labels,'Fpz')),1:2)-...
%                  obj.pts1020(find(ismember(obj.labels,'Cz')),1:2));
%             
            % add a circle for the head
            theta = linspace(0,2*pi);
            plot(headradius*cos(theta),headradius*sin(theta),'color',[.4 .4 .4],'linestyle','--');
            
            headradius=norm([x(find(ismember(obj.labels,'nas'))) y(find(ismember(obj.labels,'nas')))]-...
                [x(find(ismember(obj.labels,'Cz'))) y(find(ismember(obj.labels,'Cz')))]);
            
            plot(headradius*cos(theta),headradius*sin(theta),'k');
            plot([-headradius headradius],[0 0],'color',[.6 .6 .6],'linestyle','--')
            plot([0 0],[-headradius headradius],'color',[.6 .6 .6],'linestyle','--')
            
            line([-10 0],[-headradius -headradius-10],'color','k');
            line([10 0],[-headradius -headradius-10],'color','k');
            scatter([-15 15],[-headradius -headradius],'filled','k','sizedata',120);
            
            % Draw the central sulcus
           
            %lst={'CpZ','C3','C5-FC5'}
            pts=[x(find(ismember(obj.labels,'CPz'))), x(find(ismember(obj.labels,'C3'))), ...
                (x(find(ismember(obj.labels,'C5')))+x(find(ismember(obj.labels,'FC5'))))/2;...
                y(find(ismember(obj.labels,'CPz'))), y(find(ismember(obj.labels,'C3'))), ...
                (y(find(ismember(obj.labels,'C5')))+y(find(ismember(obj.labels,'FC5'))))/2]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            xx=[min(pts(:,1)):.1:max(pts(:,1))];
            p=polyfit(pts(:,1),pts(:,2),2);
            plot(xx,polyval(p,xx),'color',[.4 .4 .4],'linestyle','-');
            
            %lst2={'CpZ','C4','C6-FC6'}
            pts=[x(find(ismember(obj.labels,'CPz'))), x(find(ismember(obj.labels,'C4'))), ...
                (x(find(ismember(obj.labels,'C6')))+x(find(ismember(obj.labels,'FC6'))))/2;...
                y(find(ismember(obj.labels,'CPz'))), y(find(ismember(obj.labels,'C4'))), ...
                (y(find(ismember(obj.labels,'C6')))+y(find(ismember(obj.labels,'FC6'))))/2]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            xx=[min(pts(:,1)):.1:max(pts(:,1))];
            p=polyfit(pts(:,1),pts(:,2),2);
            plot(xx,polyval(p,xx),'color',[.4 .4 .4],'linestyle','-');
            
            % Add the insular sulcus
           % lst={'FT9','FT7','C5','Cp5'}
           
             pts=[x(find(ismember(obj.labels,'FT9'))), x(find(ismember(obj.labels,'FT7'))), ...
                  x(find(ismember(obj.labels,'C5'))), x(find(ismember(obj.labels,'CP5')));...
                  y(find(ismember(obj.labels,'FT9'))), y(find(ismember(obj.labels,'FT7'))), ...
                  y(find(ismember(obj.labels,'C5'))), y(find(ismember(obj.labels,'CP5')))]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            yy=[min(pts(:,2)):.1:max(pts(:,2))];
            p=polyfit(pts(:,2),pts(:,1),3);
            plot(polyval(p,yy),yy,'color',[.4 .4 .4],'linestyle','-');
            
            % lst={'FT10','FT8','C6','Cp6'}
            pts=[x(find(ismember(obj.labels,'FT10'))), x(find(ismember(obj.labels,'FT8'))), ...
                  x(find(ismember(obj.labels,'C6'))), x(find(ismember(obj.labels,'CP6')));...
                  y(find(ismember(obj.labels,'FT10'))), y(find(ismember(obj.labels,'FT8'))), ...
                  y(find(ismember(obj.labels,'C6'))), y(find(ismember(obj.labels,'CP6')))]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            yy=[min(pts(:,2)):.1:max(pts(:,2))];
            p=polyfit(pts(:,2),pts(:,1),3);
            plot(polyval(p,yy),yy,'color',[.4 .4 .4],'linestyle','-');
            
            set(gca,'YDir','reverse');
            set(gcf,'color','w');
            
            if(obj.zoom)
                dx=(max(xop)-min(xop))/10;
                dy=(max(yop)-min(yop))/10;
                set(gca,'Xlim',[min(xop)-dx max(xop)+dx]);
                set(gca,'Ylim',[-headradius*1.13 max(yop)+dy]);
            end
            set(gca,'Xdir','reverse');
            if(nargout>0)
                varargout{1}=h;
            end
            
        end
        
        function headcircum = get.headcircum(obj)
            a=obj.LR_distance/2;
            b=obj.AP_distance/2;
            headcircum=.9*pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));
            
        end
        
        function LR_arclength = get.LR_arclength(obj)
            a=obj.LR_distance/2;
            b=obj.IS_distance;
            LR_arclength=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))/2;
            
        end
        
        function AP_arclength = get.AP_arclength(obj)
            a=obj.AP_distance/2;
            b=obj.IS_distance;
            AP_arclength=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))/2;
            
        end
        
        function [x,y]=convert2d(obj,pts)
            
            %Clarke far-side general prospective azumuthal projection
            r = -2.4;
            R=sqrt(sum(pts.^2,2));
            x=r*R.*(pts(:,1)./abs(pts(:,3)-r*R));
            y=r*R.*(pts(:,2)./abs(pts(:,3)-r*R));
            
        end
        
        %draw1020image;  % Code to draw an image in 10-20space
        % makeimage;  % code to do simple image reconstruction in 10-20 space
        %registerprobe;
        % depth map
        % define ROIs
        % coregister
        % register probe
        
        
        
    end
    
end
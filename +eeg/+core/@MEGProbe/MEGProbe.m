classdef MEGProbe
    %% PROBE - This object hold EEG geometries.
    %
    % Properties:
    %    electrodes-  a table containing the position, name, and type of any
    %
    %  Methods:
    %     draw   - displays a visualization of the probe geometry
    
    properties
        electrodes    % table describing tposition of electrodes
        link
    end
    
    methods
        function obj = MEGProbe(hdr)
            %% Probe - Creates a probe object.
            %
            % Args:
            %     labels - name of 10-20 points
           
           lst=find([hdr.chs.kind]==1);
           coils={hdr.ch_names{:}};
           xyz =horzcat(hdr(:).chs.loc);
           
           xyz=xyz(1:3,:)';
           xyz(:,4)=1;
           xyz=xyz*hdr.dev_head_t.trans;
           xyz=xyz*1000;
           
           % Register to the 10-20 system
           ptsfid = [hdr.dig(1:3).r]'*1000;
           
           tbl=nirs.util.list_1020pts('?');
          
           fid=[tbl([2 1 3],:).X tbl([2 1 3],:).Y tbl([2 1 3],:).Z];
           
           cfid = mean(fid,1);
           cptsfid = mean(ptsfid,1);
           trans= cfid-cptsfid;
           fid=fid-ones(3,1)*cfid;
           ptsfid=ptsfid-ones(3,1)*cptsfid;
           
           fid(4,:)=cross(fid(1,:),fid(2,:));
           fid(5,:)=cross(fid(3,:),fid(2,:));
           ptsfid(4,:)=cross(ptsfid(1,:),ptsfid(2,:));
           ptsfid(5,:)=cross(ptsfid(3,:),ptsfid(2,:));
           
           
           for i=1:5;
               fid(i,:)=fid(i,:)/norm(fid(i,:));
               ptsfid(i,:)=ptsfid(i,:)/norm(ptsfid(i,:));
           end
           T=ptsfid\fid;
       %xyz=xyz(:,1:3)*T;
       % xyz=xyz+ones(size(xyz,1),1)*trans;
           obj.link=table;
           
           X=xyz(lst,1);
           Y=xyz(lst,2);
           Z=xyz(lst,3);
           Name = {coils{lst}}';
           for i=1:length(Name);
               if(hdr.chs(i).coil_type==3012)
                   type{i,1}='Grad';
               else
                   type{i,1}='Mag';
               end
           end
           for i=1:3:length(Name)
               type{i}=[type{i} '-A'];
           end
           for i=2:3:length(Name)
               type{i}=[type{i} '-B'];
           end
          obj.link=table(Name,type);
          obj.electrodes=table(Name,X,Y,Z);
           
                      
        end
        
        function varargout=draw(obj, colors, lineStyles, axis_handle )
            
            types= unique(obj.link.type);
            for i=1:length(types)
                lst=find(ismember(obj.link.type,types{i}));
                if nargin < 4
                    axis_handle = subplot(1,length(types),i);
                end
                
               % figure;
                
            n=height(obj.electrodes(lst,:));
            
            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [height(obj.link) 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [height(obj.link) 1]);
            end
            
            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [height(obj.link) 1]);
            elseif size(lineStyles, 1) == 1
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [height(obj.link) 1]);
            end
            
            
                     hold on;
            tbl=nirs.util.list_1020pts('?');
            pts=[tbl.X tbl.Y tbl.Z];
            r = -2.4;
            R=sqrt(sum(pts.^2,2));
            x=r*R.*(pts(:,1)./abs(pts(:,3)-r*R));
            y=r*R.*(pts(:,2)./abs(pts(:,3)-r*R));
            
            dx=-x(find(ismember(tbl.Name,'Cz')));
            dy=-y(find(ismember(tbl.Name,'Cz')));
            
            scatter(x+dx,y+dy,'filled','MarkerFaceColor',[.8 .8 .8]);
            
            pts=[obj.electrodes.X obj.electrodes.Y obj.electrodes.Z];
            R=sqrt(sum(pts.^2,2));
            xele=r*R.*(pts(:,1)./abs(pts(:,3)-r*R))+dx;
            yele=r*R.*(pts(:,2)./abs(pts(:,3)-r*R))+dy;
            
            for i=1:height(obj.electrodes(lst,:))
                %h(iChan) = line(x, y, 'Color', colors(iChan, :), lineStyles{iChan, :});
                if(strcmp(lineStyles{lst(i), 2},'--'))
                    h(lst(i))=scatter(xele(lst(i)),yele(lst(i)),'filled','MarkerFaceColor',colors(lst(i), :));
                    set(h(lst(i)),'SizeData',50);
                else
                    h(lst(i))=scatter(xele(lst(i)),yele(lst(i)),'filled','MarkerFaceColor',colors(lst(i), :),'MarkerEdgeColor','k');
                    set(h(lst(i)),'SizeData',100);
                end
                t(lst(i))=text(xele(lst(i)),yele(lst(i)),obj.electrodes.Name{lst(i)});
            end
            set(t(lst),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
            
          headradius=norm([x(find(ismember(tbl.Name,'Fpz'))) y(find(ismember(tbl.Name,'Fpz')))]-...
                [x(find(ismember(tbl.Name,'Cz'))) y(find(ismember(tbl.Name,'Cz')))]);
            
            theta = linspace(0,2*pi);
            plot(headradius*cos(theta),headradius*sin(theta),'color',[.4 .4 .4],'linestyle','--');
            
            headradius=norm([x(find(ismember(tbl.Name,'nas'))) y(find(ismember(tbl.Name,'nas')))]-...
                [x(find(ismember(tbl.Name,'Cz'))) y(find(ismember(tbl.Name,'Cz')))]);
            
            plot(headradius*cos(theta),headradius*sin(theta),'k');
            plot([-headradius headradius],[0 0],'color',[.6 .6 .6],'linestyle','--')
            plot([0 0],[-headradius headradius],'color',[.6 .6 .6],'linestyle','--')
            
            line([-10 0],[-headradius -headradius-10],'color','k');
            line([10 0],[-headradius -headradius-10],'color','k');
            scatter([-15 15],[-headradius -headradius],'filled','k','sizedata',120);
            
            % Draw the central sulcus
            
            %lst={'CpZ','C3','C5-FC5'}
            pts=[x(find(ismember(tbl.Name,'CPz'))), x(find(ismember(tbl.Name,'C3'))), ...
                (x(find(ismember(tbl.Name,'C5')))+x(find(ismember(tbl.Name,'FC5'))))/2;...
                y(find(ismember(tbl.Name,'CPz'))), y(find(ismember(tbl.Name,'C3'))), ...
                (y(find(ismember(tbl.Name,'C5')))+y(find(ismember(tbl.Name,'FC5'))))/2]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            xx=[min(pts(:,1)):.1:max(pts(:,1))];
            p=polyfit(pts(:,1),pts(:,2),2);
            plot(xx,polyval(p,xx),'color',[.4 .4 .4],'linestyle','-');
            
            %lst2={'CpZ','C4','C6-FC6'}
            pts=[x(find(ismember(tbl.Name,'CPz'))), x(find(ismember(tbl.Name,'C4'))), ...
                (x(find(ismember(tbl.Name,'C6')))+x(find(ismember(tbl.Name,'FC6'))))/2;...
                y(find(ismember(tbl.Name,'CPz'))), y(find(ismember(tbl.Name,'C4'))), ...
                (y(find(ismember(tbl.Name,'C6')))+y(find(ismember(tbl.Name,'FC6'))))/2]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            xx=[min(pts(:,1)):.1:max(pts(:,1))];
            p=polyfit(pts(:,1),pts(:,2),2);
            plot(xx,polyval(p,xx),'color',[.4 .4 .4],'linestyle','-');
            
            % Add the insular sulcus
            % lst={'FT9','FT7','C5','Cp5'}
            
            pts=[x(find(ismember(tbl.Name,'FT9'))), x(find(ismember(tbl.Name,'FT7'))), ...
                x(find(ismember(tbl.Name,'C5'))), x(find(ismember(tbl.Name,'CP5')));...
                y(find(ismember(tbl.Name,'FT9'))), y(find(ismember(tbl.Name,'FT7'))), ...
                y(find(ismember(tbl.Name,'C5'))), y(find(ismember(tbl.Name,'CP5')))]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            yy=[min(pts(:,2)):.1:max(pts(:,2))];
            p=polyfit(pts(:,2),pts(:,1),3);
            plot(polyval(p,yy),yy,'color',[.4 .4 .4],'linestyle','-');
            
            % lst={'FT10','FT8','C6','Cp6'}
            pts=[x(find(ismember(tbl.Name,'FT10'))), x(find(ismember(tbl.Name,'FT8'))), ...
                x(find(ismember(tbl.Name,'C6'))), x(find(ismember(tbl.Name,'CP6')));...
                y(find(ismember(tbl.Name,'FT10'))), y(find(ismember(tbl.Name,'FT8'))), ...
                y(find(ismember(tbl.Name,'C6'))), y(find(ismember(tbl.Name,'CP6')))]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            yy=[min(pts(:,2)):.1:max(pts(:,2))];
            p=polyfit(pts(:,2),pts(:,1),3);
            plot(polyval(p,yy),yy,'color',[.4 .4 .4],'linestyle','-');
            
            set(gca,'YDir','reverse');
            set(gcf,'color','w');
            set(gca,'Xdir','reverse');
            axis off;
            end
            if(nargout==1)
                varargout{1}=h;
            end
        end
    end
end


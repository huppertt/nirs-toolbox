classdef Probe
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
        function obj = Probe( labels)
            %% Probe - Creates a probe object.
            %
            % Args:
            %     labels - name of 10-20 points
            
            tbl=nirs.util.list_1020pts('?');
            lst=find(ismember(lower(tbl.Name),lower(labels)));
            tbl.Type=[];
           
            obj.electrodes=tbl(lst,:);
            
            electrode=[1:length(lst)]';
            type=repmat({'eeg'},length(lst),1);
            obj.link=table(electrode,type);
            
        end
        
        function varargout=draw(obj, colors, lineStyles, axis_handle )
            
            n=height(obj.electrodes);
            
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
            
            for i=1:height(obj.electrodes)
                %h(iChan) = line(x, y, 'Color', colors(iChan, :), lineStyles{iChan, :});
                if(strcmp(lineStyles{i, 2},'--'))
                    h(i)=scatter(xele(i),yele(i),'filled','MarkerFaceColor',colors(i, :));
                    set(h(i),'SizeData',50);
                else
                    h(i)=scatter(xele(i),yele(i),'filled','MarkerFaceColor',colors(i, :),'MarkerEdgeColor','k');
                    set(h(i),'SizeData',100);
                end
                t(i)=text(xele(i),yele(i),obj.electrodes.Name{i});
            end
            set(t,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12);
            
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
            if(nargout==1)
                varargout{1}=h;
            end
        end
    end
end


function varargout=draw( obj, colors, lineStyles, axis_handle)
    %% draw - Plots the probe geometry.
    % 
    % Args:
    %     colors      - (optional) n x 3 array of colors [R, G, B] for each channel
    %     lineStyles  - (optional) 2D cell array with each row containing
    %                   arguments for the 'line' functions (e.g. {'LineWidth',6})
    %     axis_handle - (optional) handle to axis to the plot to
        
    % sd pairs
    if nargin < 2
        colors=[];
    end
    if nargin < 3 
        lineStyles = [];
    end
    if nargin < 4
        axis_handle = axes();
    end
    
    p=obj.convertProbe1020;
    hp=p.draw(colors, lineStyles, axis_handle);

    link_all=obj.link;
    link_all=link_all(ismember(link_all.type,link_all.type(1)),:);
    SrcA=cellfun(@(x)str2num(x(3:strfind(x,':')-1)),link_all.source);
    SrcB=cellfun(@(x)str2num(x(3:strfind(x,':')-1)),link_all.detector);
    DetA=cellfun(@(x)str2num(x(strfind(x,'D-')+2:end)),link_all.source);
    DetB=cellfun(@(x)str2num(x(strfind(x,'D-')+2:end)),link_all.detector);

    link=obj.link_probe;
    
    xyz=zeros(length(hp),3);
    for i=1:length(hp)
        xyz(i,1)=mean(get(hp(i),'Xdata'));
        xyz(i,2)=mean(get(hp(i),'Ydata'));
        if(isempty(get(hp(i),'ZData')))
            xyz(i,3)=NaN;
        else
            xyz(i,3)=mean(get(hp(i),'Zdata'));
        end

    end
    tbl=p.link;
    tbl=tbl(ismember(tbl.type,tbl.type(1)),:);

    for i=1:size(SrcA,1) 
        ee(i,1)=find(SrcB(i)==tbl.source &...
                    DetB(i)==tbl.detector); 
        ss(i,1)=find(SrcA(i)==tbl.source &...
                    DetA(i)==tbl.detector); 
        
    end;
  
    h = drawConnections(xyz(ss,:),xyz(ee,:),[SrcA DetA SrcB DetB],axis_handle);
    set(h,'color',[.7 .7 .7]);
    
    if(nargout==1)
        varargout{1}=h;
    end
    
end

function h=drawConnections(st, en,ll, axis_handle)
    hold(axis_handle,'on');
    for iChan=1:size(st,1)
        if(isnan(st(iChan,3)))
            h(iChan)=line([st(iChan,1) en(iChan,1)],[st(iChan,2) en(iChan,2)]);
        else
            h(iChan)=line([st(iChan,1) en(iChan,1)],[st(iChan,2) en(iChan,2)],[st(iChan,3) en(iChan,3)]);
        end
        set(h(iChan),'UserData',ll(iChan,:));
    end
    hold(axis_handle,'off');
end

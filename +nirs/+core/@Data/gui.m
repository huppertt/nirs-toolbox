function h = gui(obj)
% this launchs the NIRSviewer GUI for this data

h=figure('Tag','nirsviewer');
set(h,'MenuBar','none');

handles=guihandles(h);
handles.axis_main = axes('parent',h,...
    'Units','normalized','Box','on','Position',[.05 .1 .6 .8]);
xlabel(handles.axis_main,'time (s)');

handles.axis_sdg = axes('parent',h,...
    'Units','normalized','Box','on','Position',[.7 .55 .28 .35],...
    'XTick',[],'YTick',[]);

try
    types=arrayfun(@(x){num2str(x)},obj(1).probe.link.type);
catch
    types=arrayfun(@(x){x{1}},obj(1).probe.link.type);
end

handles.selecttype = uicontrol('Style','popupmenu','String',unique(types),...
    'Units','normalized','Position',[.69 .43 .31 .1]);
set(handles.selecttype,'callback',@changetype);

linehandles=obj(1).draw(1:height(obj(1).probe.link),[],handles.axis_main);
set(linehandles,'tag','dataline');

SDcolors=nirs.util.makeSDcolors(obj(1).probe.link);

set(linehandles,'visible','off');

for i=1:length(linehandles)
    set(linehandles(i),'userdata',types{i});
end


subtype=get(handles.selecttype,'String');

for idx=1:length(subtype)
    lst=find(ismember(types,subtype{idx}));
    set(linehandles(lst),'visible','off');
    for idx2=1:length(lst)
        set(linehandles(lst(idx2)),'color',SDcolors(idx2,:));
    end
    LstAll(:,idx)=lst;
end

subtype={subtype{get(handles.selecttype,'value')}};

for idx=1:length(subtype)
    lst=find(ismember(types,subtype{idx}));
    set(linehandles(lst),'visible','on');
    for idx2=1:length(lst)
        set(linehandles(lst(idx2)),'color',SDcolors(idx2,:));
    end
end

SDGhandlesBase=obj(1).probe.draw([],[],handles.axis_sdg);
set(SDGhandlesBase,'Color',[.8 .8 .8]);

lines=findobj('type','line','parent',handles.axis_main,'tag','dataline');
isvis={};
for idx=1:length(lines)
    isvis{idx}=get(lines(idx),'visible');
end

SDGhandles=obj(1).probe.draw([],[],handles.axis_sdg);
ii=get(handles.selecttype,'value');
for idx=1:length(SDGhandles)
    set(SDGhandles(idx),'tag',['SDG' num2str(idx)]);
    set(SDGhandlesBase(idx),'tag',['SDG' num2str(idx)]);
    set(SDGhandles(idx),'color',SDcolors(idx,:));
    linelinks(idx)=linkprop([linehandles(LstAll(idx,ii)),SDGhandles(idx)],{'Visible','Color'});
end


set(SDGhandles,'ButtonDownFcn','set(gcbo,''visible'',''off'');');
set(SDGhandlesBase,'ButtonDownFcn','set(findobj(''tag'',get(gcbo,''tag'')),''visible'',''on'');');
setappdata(h,'linelinks',linelinks);
axis(handles.axis_sdg,'tight');

handles.SDG_menu=uicontextmenu;
uimenu(handles.SDG_menu,'Label','All Off','Callback',@setalloff);
uimenu(handles.SDG_menu,'Label','All On','Callback',@setallon);

set(handles.axis_main,'uiContextMenu',handles.SDG_menu);

if(isempty(nirs.getStimNames(obj)))
    legend off;
end


handles.LstAll=LstAll;
handles.linelinks=linelinks;
handles.SDGhandles=SDGhandles;
handles.linehandles=linehandles;

if(length(obj)>1)
    warning('off','MATLAB:Java:DuplicateClass');
    handles.uipanel=uipanel('parent',h,...
        'Units','normalized','Position',[.7 .05 .28 .4]);
    
    demo=nirs.createDemographicsTable(obj);
    if(isempty(demo))
        DataTable=cell(length(obj),1);
        for i=1:length(obj)
            if(isempty(obj(i).description))
                DataTable{i}=['Entry-' num2str(i)];
            else
                [~, DataTable{i},~]=fileparts(obj(i).description);
            end
        end
        jtable = treeTable(handles.uipanel,{'Entry'},DataTable,...
            'ColumnTypes',{'char'},'groups',[false]);
    else
        cnt=length(find(ismember({'group','subject'},demo.Properties.VariableNames)));
        DataTable=cell(length(obj),1+cnt);
        for i=1:height(demo)
            if(isempty(obj(i).description))
                DataTable{i,1}=['Entry-' num2str(i)];
            else
                [~, DataTable{i,1},~]=fileparts(obj(i).description);
            end
            if(cnt==2)
                DataTable{i,2}=demo.group{i};
                DataTable{i,3}=demo.subject{i};
                str={'Entry' 'Group' 'Subject'};
            elseif(ismember({'group'},demo.Properties.VariableNames))
                DataTable{i,2}=demo.group{i};
                str={'Entry' 'Group'};
            elseif(ismember({'subject'},demo.Properties.VariableNames))
                DataTable{i,2}=demo.subject{i};
                
                str={'Entry'  'Subject'};
            else
                str={'Entry'};
            end
        end
        jtable = treeTable(handles.uipanel,str,DataTable,...
            'ColumnTypes',repmat({'char'},cnt+1,1),'groups',false(cnt+1,1));
    end
    
    set(handle(jtable.getSelectionModel,'CallbackProperties'), 'ValueChangedCallback', []);
    set(jtable,'MousePressedCallback',@tablechange);
end

set(handles.axis_main,'userdata',obj);

set(h,'userdata',handles);


return

function tablechange(varargin)
idx=get(varargin{1},'SelectedRow')+1;


handles=guihandles(findobj('tag','nirsviewer'));
handles=get(handles(1).nirsviewer,'userdata');
data=get(handles.axis_main,'userdata');

flag=false(size(handles.linehandles));
for i=1:length(handles.linehandles)
    flag(i)=strcmp(get(handles.linehandles(i),'visible'),'on');
    colors{i}=get(handles.linehandles(i),'color');
end
delete(handles.linehandles);
cla(handles.axis_main);
handles.linehandles=data(idx).draw(1:height(data(idx).probe.link),[],handles.axis_main);
for i=1:length(handles.linehandles)
    set(handles.linehandles(i),'color',colors{i});
end


delete(handles.linelinks);

set(handles.linehandles,'visible','off');

ii=get(handles.selecttype,'value');
for idx=1:length(handles.SDGhandles)
    if(any(flag(handles.LstAll(idx,:))))
        set(handles.SDGhandles(idx),'visible','on');
        set(handles.linehandles(handles.LstAll(idx,ii)),'visible','on');
    end
    
    handles.linelinks(idx)=linkprop([handles.linehandles(handles.LstAll(idx,ii)),handles.SDGhandles(idx)],{'Visible','Color'});
end
set(handles.nirsviewer,'userdata',handles);


return



function setalloff(varargin)
handles=guihandles(gcbo);
handles=get(handles.nirsviewer,'userdata');
 set(handles.SDGhandles,'visible','off');
return

function setallon(varargin)
handles=guihandles(gcbo);
handles=get(handles.nirsviewer,'userdata');
 set(handles.SDGhandles,'visible','on');
return



function changetype(varargin)

handles=guihandles(gcbo);
handles=get(handles.nirsviewer,'userdata');
subtype=get(handles.selecttype);

flag=false(length(handles.linehandles),1);
for i=1:length(handles.linehandles)
   flag(i)=strcmp(get(handles.linehandles(i),'visible'),'on');
end

delete(handles.linelinks);

set(handles.linehandles,'visible','off');

ii=get(handles.selecttype,'value');
for idx=1:length(handles.SDGhandles)
    if(any(flag(handles.LstAll(idx,:))))
        set(handles.SDGhandles(idx),'visible','on');
        set(handles.linehandles(handles.LstAll(idx,ii)),'visible','on');
    end
    
    handles.linelinks(idx)=linkprop([handles.linehandles(handles.LstAll(idx,ii)),handles.SDGhandles(idx)],{'Visible','Color'});
end
set(handles.nirsviewer,'userdata',handles);

return

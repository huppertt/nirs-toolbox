function varargout = nirs_viewer(varargin)
% NIRS_VIEWER MATLAB code for nirs_viewer.fig

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nirs_viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @nirs_viewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before nirs_viewer is made visible.
function nirs_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);



return

% --- Outputs from this function are returned to the command line.
function varargout = nirs_viewer_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
return

function varargout = nirs_viewer_dataview(varargin) 

handles=guihandles(findobj('tag','figure_nirsview'));
name=get(handles.listbox_data,'String');
try
    name=name{get(handles.listbox_data,'value')};
catch
    %Nothing loaded yet
    return
end
subtype=strtrim(name(strfind(name,':')+1:end));
name=strtrim(name(1:strfind(name,':')-1));
subtypesAll=evalin('base',['unique(' name '(1).probe.link.type);']);

if(~iscell(subtypesAll)); subtypesAll=num2cell(subtypesAll); end;
for idx2=1:length(subtypesAll)
    if(isnumeric(subtypesAll{idx2})); subtypesAll{idx2}=num2str(subtypesAll{idx2}); end;
end;
set(handles.listbox_showdatasubtype,'String',subtypesAll,'value',find(ismember(subtype,subtypesAll)));

data=evalin('base',name);
handles=guihandles(findobj('tag','figure_nirsview'));


%% Populate the filename/subject/group info

table=nirs.createDemographicsTable(data);
if(isfield(table,'Group')); 
    table.group=table.Group;
    table=rmfield(table,'Group');
end
if(isfield(table,'Subject')); 
    table.subject=table.Subjects;
    table=rmfield(table,'Subject');
end

if(~isfield(table,'group'))
    for idx=1:length(data)
        table.group{idx}='group 1';
    end
end
if(~isfield(table,'subject'))
    for idx=1:length(data)
        table.subject{idx}='subject 1';
    end
end

delete(findobj('tag','uitree_cont'));
delete(findobj('tag','uitree_obj'));

import javax.swing.*
import javax.swing.tree.*;
root = uitreenode('v0','Subjects', 'Subjects', [], false);

groups=unique(table.group);
for gIdx=1:length(groups)
    g = uitreenode('v0', groups{gIdx},  groups{gIdx}, [], false);
    lst=find(ismember(table.group, groups{gIdx}));
    subj=table.subject;
    subj=unique({subj{lst}});
    for sIdx=1:length(subj)

        s = uitreenode('v0', subj{sIdx},  subj{sIdx}, [], false);
        lstfiles=find(ismember(table.group, groups{gIdx}) & ismember(table.subject, subj{sIdx}));
        for cf=1:length(lstfiles)
            [~,file,ext]=fileparts(data(lstfiles(cf)).description);
            f = uitreenode('v0', num2str(lstfiles(cf)),[file ext],[],true);
            s.add(f);
        end
        g.add(s);
    end
    root.add(g);
end

[mtree,container] = uitree('v0', 'Root', root);
set(container,'units',get(handles.uipanel1,'units'));
set(container,'tag','uitree_cont');

set(container,'position',get(handles.uipanel1,'position'));
set(container,'visible','on');
set(container,'userdata',mtree);

mtree.expand(mtree.getRoot);
mtree.expand(s);
mtree.setSelectedNode(f);
set(mtree,'NodeSelectedCallback',@updatewindow);

updatewindow;
return

% This function does the drawing 
function updatewindow(varargin)

handles=guihandles(findobj('tag','figure_nirsview'));
name=get(handles.listbox_data,'String');
name=name{get(handles.listbox_data,'value')};
subtype={strtrim(name(strfind(name,':')+1:end))};
name=strtrim(name(1:strfind(name,':')-1));
data=evalin('base',name);


a=findobj('tag','uitree_cont');
a=get(a,'Userdata');
node=get(a.Tree,'LastSelectedPathComponent');

filename=node.getName;
val=str2num(node.getValue);
if(isempty(val))
    return
end

try
    typesAll=arrayfun(@(x){num2str(x)},data(val).probe.link.type);
catch
    typesAll=arrayfun(@(x){x{1}},data(val).probe.link.type);
end

cla(handles.axes_SDG);
axes(handles.axes_SDG);
hold on;
SDGhandlesBase=data(val).probe.draw([],[],handles.axes_SDG);
set(SDGhandlesBase,'Color',[.8 .8 .8]);

SDGhandles=data(val).probe.draw([],[],handles.axes_SDG);
axis tight;
axis on;
box off;
set(handles.axes_SDG,'color',get(handles.figure_nirsview,'color'))
set(handles.axes_SDG,'Xcolor',get(handles.figure_nirsview,'color'))
set(handles.axes_SDG,'Ycolor',get(handles.figure_nirsview,'color'))

set(handles.axes_SDG,'Xtick',[],'Ytick',[])

%See if there is anything already drawn and keep the visibility (if
%possible)
lines=findobj('type','line','parent',handles.axes_maindata,'tag','dataline');
isvis={};
for idx=1:length(lines)
    isvis{idx}=get(lines(idx),'visible');
end

cla(handles.axes_maindata);
axes(handles.axes_maindata);
linehandles=data(val).draw;
set(linehandles,'tag','dataline');

SDcolors=nirs.util.makeSDcolors(data(val).probe.link);

set(linehandles,'visible','off');
for idx=1:length(subtype)
    lst=find(ismember(typesAll,subtype{idx}));
    set(linehandles(lst),'visible','on');
    for idx2=1:length(lst)
        set(linehandles(lst(idx2)),'color',SDcolors(idx2,:));
    end
    LstAll(:,idx)=lst;
end

for idx=1:length(SDGhandles)
    set(SDGhandles(idx),'tag',['SDG' num2str(idx)]);
    set(SDGhandlesBase(idx),'tag',['SDG' num2str(idx)]);
    set(SDGhandles(idx),'color',SDcolors(idx,:));
    linelinks(idx)=linkprop([linehandles(LstAll(idx,:)),SDGhandles(idx)],{'Visible','Color'});
end
for idx=1:size(LstAll,2)
    linelinks(end+1)=linkprop(linehandles(LstAll(:,idx)),'LineStyle');
end

set(SDGhandles,'ButtonDownFcn','set(gcbo,''visible'',''off''); nirs_viewer(''updatewin'');');
set(SDGhandlesBase,'ButtonDownFcn','set(findobj(''tag'',get(gcbo,''tag'')),''visible'',''on''); nirs_viewer(''updatewin'');');
setappdata(handles.figure_nirsview,'linelinks',linelinks);

set(handles.axes_maindata,'uiContextMenu',handles.MainPlot_menu);
set(handles.axes_SDG,'uiContextMenu',handles.SDG_menu);

delete(linehandles(find(~ismember(typesAll,subtype))));

lines=findobj('type','line','parent',handles.axes_maindata,'tag','dataline');
if(length(lines)==length(isvis))
    for idx=1:length(lines)
        set(lines(idx),'visible',isvis{idx});
    end
end

return

function updatedatalist
handles=guihandles(findobj('tag','figure_nirsview'));

% Find all the nirs.core.data types avaliable for display
datatypes=evalin('base','whos;');
datatypes(~ismember({datatypes.class},'nirs.core.Data'))=[];
names={datatypes.name};
subnames={};
for idx=1:length(names)
    subtypes=evalin('base',['unique(' names{idx} '(1).probe.link.type);']);
    if(~iscell(subtypes)); subtypes=num2cell(subtypes); end;
    for idx2=1:length(subtypes)
        if(isnumeric(subtypes{idx2})); subtypes{idx2}=num2str(subtypes{idx2}); end;
        subnames={subnames{:} [names{idx} ' : ' subtypes{idx2}]};
    end
end

if(length(datatypes)>0)
    set(handles.listbox_data,'String',subnames);
else
    set(handles.listbox_data,'String','<None Loaded>');  
end


%Update the UImenu with the list of avaliable channel stats objected
datatypes=evalin('base','whos;');
datatypes(~ismember({datatypes.class},'nirs.core.ChannelStats'))=[];
names={datatypes.name};
delete(get(handles.uimenu_launch_chanstats,'children'));
set(handles.uimenu_launch_chanstats,'Callback',[]);
for idx=1:length(datatypes)
    h=uimenu('parent',handles.uimenu_launch_chanstats,'Label',datatypes(idx).name);
    set(h,'Callback',['nirs_viewer(''uimenu_launch_chanstats_Callback'',gcbo,[],guidata(gcbo),''' datatypes(idx).name ''')']);
end

%% TODO-  the HTML reporter needs to handle core.data class
if(0)
datatypes2=evalin('base','whos;');
datatypes2(~ismember({datatypes2.class},'nirs.core.Data'))=[];
names={datatypes.name datatypes2.name};
delete(get(handles.uimenu_create_datareport,'children'));
set(handles.uimenu_create_datareport,'Callback',[]);
for idx=1:length(names)
    h=uimenu('parent',handles.uimenu_create_datareport,'Label',names{idx});
    set(h,'Callback',['nirs_viewer(''uimenu_create_datareport_Callback'',gcbo,[],guidata(gcbo),''' names{idx} ''')']);
end
end

return

function updatewin
updatedatalist;
 
handles=guihandles(findobj('tag','figure_nirsview'));

if(strcmp(get(handles.uimenu_autoscale,'Checked'),'on'))
    l=findobj('type','line','parent',handles.axes_maindata,'visible','on','Tag','dataline');
    lstim=findobj('type','line','parent',handles.axes_maindata,'visible','on','Tag','');
else
    l=findobj('type','line','parent',handles.axes_maindata,'Tag','dataline');
    lstim=findobj('type','line','parent',handles.axes_maindata,'Tag','');
end

rangey=[-1 1];


for idx=1:length(l)
    rangey(idx,1)=min(get(l(idx),'Ydata'));
    rangey(idx,2)=max(get(l(idx),'Ydata'));
    rangex(idx,1)=min(get(l(idx),'Xdata'));
    rangex(idx,2)=max(get(l(idx),'Xdata'));
end

minY=min(rangey(:));
rangeY = max(rangey(:))-minY;

for idx=1:length(lstim)
    yd=get(lstim(idx),'yData');
    yd=yd-min(yd);
    yd=(yd./max(yd))*rangeY*.1+minY-rangeY*.15;
    set(lstim(idx),'yData',yd);
    rangex(length(l)+idx,1)=min(get(lstim(idx),'Xdata'));
   rangex(length(l)+idx,2)=max(get(lstim(idx),'Xdata'));
end
minY=minY-rangeY*.15;


set(handles.axes_maindata,'Ylim',[minY max(rangey(:))],'Xlim',[min(rangex(:)) max(rangex(:))]);

return

% --- Executes on selection change in listbox_data.
function listbox_data_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_data
name=get(hObject,'String');

if(strcmp(name,'<None Loaded>'))
    return
end

name=name{get(hObject,'value')};
nirs_viewer_dataview(name);

% --- Executes during object creation, after setting all properties.
function listbox_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
updatedatalist;
nirs_viewer_dataview;

return

% --------------------------------------------------------------------
function uimenu_about_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uimenu_helpwin_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_helpwin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('https://bitbucket.org/jeffx/nirs-toolbox/wiki/Home','-browser')


% --------------------------------------------------------------------
function uimenu_jobmanager_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_jobmanager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=nirs.viz.jobmanager;

set(h,'CloseRequestFcn',[get(h,'CloseRequestFcn') '; nirs_viewer(''uimenu_refresh_Callback'');']);

return

% --------------------------------------------------------------------
function uimenu_stimdesign_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_stimdesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uimenu_editdemos_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_editdemos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uimenu_loadfiles_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_loadfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=nirs.io.loadGUI;
waitfor(h);
uimenu_refresh_Callback([],[],[]);

% --------------------------------------------------------------------
function uimenu_loadsaved_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_loadsaved (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Load a saved workspace');
if(filename~=0)
    evalin('base',['load(''' fullfile(pathname,filename) ''');']);
end
uimenu_refresh_Callback([],[],[]);


return

% --------------------------------------------------------------------
function uimenu_savedata_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Find all the variables associated with the nirs.core toolbox and save
% those to file
datatypes=evalin('base','whos;');
names={};
for idx=1:length(datatypes)
    if(~isempty(strfind(datatypes(idx).class,'nirs.')))
        names={names{:} datatypes(idx).name};
    end
end

if(~isempty(names))
    [filename, pathname] = uiputfile('nirs_results.mat', 'Save Workspace as');
    if(filename~=0)
        filename=fullfile(pathname,filename);
        str=['save(''' filename ''''];
        for idx=1:length(names)
            str=[str ',''' names{idx} ''''];
        end
        str=[str ');'];
        evalin('base',str);
    end
    
else
    msgbox('Nothing to save');
end


return


% --------------------------------------------------------------------
function uimenu_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatedatalist;
nirs_viewer_dataview;


% --- Executes on selection change in listbox_showdatasubtype.
function listbox_showdatasubtype_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_showdatasubtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_showdatasubtype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_showdatasubtype


% --- Executes during object creation, after setting all properties.
function listbox_showdatasubtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_showdatasubtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uimenu_autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(handles.uimenu_autoscale,'Checked'),'on'))
    set(handles.uimenu_autoscale,'Checked','off');
else
    set(handles.uimenu_autoscale,'Checked','on');
end
updatewin;
return

% --------------------------------------------------------------------
function uimenu_plotnew_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_plotnew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Launch a new figure and copy all the objects over 
f=figure;
l=legend(handles.axes_maindata);
c=copyobj([handles.axes_maindata l],f);
copyobj(handles.axes_SDG,f);

return

% --------------------------------------------------------------------
function uimenu_showall_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_showall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lines=findobj('type','line','parent',handles.axes_maindata,'tag','dataline');
set(lines,'visible','on')


return

% --------------------------------------------------------------------
function uimenu_shownone_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_shownone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lines=findobj('type','line','parent',handles.axes_maindata,'tag','dataline');
set(lines,'visible','off')

return
% --------------------------------------------------------------------
function uimenu_launch_chanstats_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to uimenu_launch_chanstats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uimenu_create_datareport_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to uimenu_create_datareport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data=evalin('base',varargin{1});
if(length(data)>1)
     ButtonName = questdlg(['Warning: This Channel Stats variable has ' num2str(length(data)) ' entries. Continue?'], ...
                         'This could take a while', ...
                         'Yes','No','No');
                     if(strcmp(ButtonName,'No'))
                         return;
                     end
end
rpt=nirs.util.create_chanstats_rpt(data,[varargin{1} '_report']);
report(rpt);

return

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

data=evalin('base',varargin{1});
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


set(handles.figure_nirsview,'Userdata',data);
updatewindow;
return

function updatewindow(varargin)


handles=guihandles(findobj('tag','figure_nirsview'));

a=findobj('tag','uitree_cont');
a=get(a,'Userdata');
node=get(a.Tree,'LastSelectedPathComponent');

filename=node.getName;
val=str2num(node.getValue);
data=get(handles.figure_nirsview,'Userdata');

cla(handles.axes_SDG);
axes(handles.axes_SDG);
hold on;
SDGhandlesBase=data(val).probe.draw([],[],handles.axes_SDG);
set(SDGhandlesBase,'Color',[.8 .8 .8]);

SDGhandles=data(val).probe.draw([],[],handles.axes_SDG);
axis tight

cla(handles.axes_maindata);
axes(handles.axes_maindata);
linehandles=data(val).draw;

typesAll=arrayfun(@(x){num2str(x)},data(val).probe.link.type);
types=unique(typesAll);
lstdisp=[1];

SDcolors=nirs.util.makeSDcolors(data(val).probe.link);

set(linehandles,'visible','off');
for idx=1:length(lstdisp)
    lst=find(ismember(typesAll,types{lstdisp(idx)}));
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
setappdata(handles.figure_nirsview,'linelinks',linelinks)


function updatewin
handles=guihandles(findobj('tag','figure_nirsview'));
datatypes=evalin('base','whos;');
datatypes(~ismember({datatypes.class},'nirs.core.Data'))=[];
set(handles.listbox_data,'String',{datatypes.name});

if(1)
    l=findobj('type','line','parent',handles.axes_maindata,'visible','on');
    for idx=1:length(l)
        rangey(idx,1)=min(get(l(idx),'Ydata'));
        rangey(idx,2)=max(get(l(idx),'Ydata'));
        rangex(idx,1)=min(get(l(idx),'Xdata'));
        rangex(idx,2)=max(get(l(idx),'Xdata'));
    end
    set(handles.axes_maindata,'Ylim',[min(rangey(:)) max(rangey(:))],'Xlim',[min(rangex(:)) max(rangex(:))]);
    
end
return

% --- Executes on selection change in listbox_data.
function listbox_data_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_data
name=get(hObject,'String');
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

% Find all the nirs.core.data types avaliable for display
datatypes=evalin('base','whos;');
datatypes(~ismember({datatypes.class},'nirs.core.Data'))=[];

if(length(datatypes)>0)
    set(hObject,'String',{datatypes.name});
    nirs_viewer_dataview(datatypes(get(hObject,'value')).name);
end

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


% --------------------------------------------------------------------
function uimenu_loadsaved_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_loadsaved (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uimenu_savedata_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

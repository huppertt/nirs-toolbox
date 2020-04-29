function varargout = StimUtil_GUI(varargin)
% STIMUTIL_GUI MATLAB code for StimUtil_GUI.fig
%      STIMUTIL_GUI, by itself, creates a new STIMUTIL_GUI or raises the existing
%      singleton*.
%
%      H = STIMUTIL_GUI returns the handle to a new STIMUTIL_GUI or the handle to
%      the existing singleton*.
%
%      STIMUTIL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STIMUTIL_GUI.M with the given input arguments.
%
%      STIMUTIL_GUI('Property','Value',...) creates a new STIMUTIL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StimUtil_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StimUtil_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StimUtil_GUI

% Last Modified by GUIDE v2.5 28-Aug-2019 12:37:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StimUtil_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @StimUtil_GUI_OutputFcn, ...
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


% --- Executes just before StimUtil_GUI is made visible.
function StimUtil_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StimUtil_GUI (see VARARGIN)

% Choose default command line output for StimUtil_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes StimUtil_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if(nargin>3)
    raw=varargin{1};
else
    raw=nirs.core.Data;
    raw(:)=[];
end

if(nargin>4)
    showdata=varargin{2};
    if(showdata)
        handles.showfNIRS.Checked='on';
    else
        handles.showfNIRS.Checked='off';
    end
else
    handles.showfNIRS.Checked='off';
end


Filenames={};
for i=1:length(raw)
    if(~isempty(raw(i).description) && exist(raw(i).description,'file')==2)
        [~,Filenames{i}]=fileparts(raw(i).description);
    elseif(~isempty(raw(i).description) && exist(raw(i).description,'file')==7)
        % NIRx data
        lst=strfind(raw(i).description,filesep);
        Filenames{i}=raw(i).description(lst(end-2)+1:lst(end)-1);
    else
        
        Filenames{i}=['File-' num2str(i)];
        raw(i).description=Filenames{i};
    end
end
set(handles.popupmenu1,'String',Filenames);
set(handles.popupmenu1,'Value',1);
set(handles.figure1,'UserData',raw);

StimUtil_GUI('updateDraw');
waitfor(handles.popupmenu1);
return

function updateDraw(varargin)

handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');

stimnames=nirs.getStimNames(raw(selected));
if(~isempty(stimnames))
if(get(handles.listbox1,'Value')>length(stimnames))
    set(handles.listbox1,'Value',1);
end

set(handles.listbox1,'String',stimnames);
idx=get(handles.listbox1,'Value');

stimevent = raw(selected).stimulus(stimnames{idx});
C={};
for i=1:length(stimevent.onset)
    C{i,1}=stimevent.name;
    C{i,2}=stimevent.onset(i);
    C{i,3}=stimevent.amp(i);
    C{i,4}=num2str(stimevent.dur(i));
end
set(handles.uitable1,'Data',C);
else
    set(handles.uitable1,'Data',{});
end

axes(handles.axes1);
cla(handles.axes1);

if(strcmp(handles.showStimMarks.Checked,'on') & ~strcmp(handles.showfNIRS.Checked,'on'))
    raw(selected).draw([]);
elseif(strcmp(handles.showStimMarks.Checked,'on') & strcmp(handles.showfNIRS.Checked,'on'))
    raw(selected).draw;
elseif(~strcmp(handles.showStimMarks.Checked,'on') & strcmp(handles.showfNIRS.Checked,'on'))
    plot(raw(selected).time,raw(selected).data);
end





return

% --- Outputs from this function are returned to the command line.
function varargout = StimUtil_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

if(~isempty(handles))
    varargout{1} = get(handles.figure1,'UserData');
else
    %canceled
    varargout{1} =[];
end
try;
    delete(handles.figure1);
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
updateDraw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
updateDraw(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uimenu_setdur_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_setdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');
stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');
stimevent = raw(selected).stimulus(stimnames{idx});

ANSWER = inputdlg('Enter duration','Duration',1,{'1'});
if(~isempty(ANSWER) && ~isempty(str2num(ANSWER{1})))
    value=str2num(ANSWER{1});
    
    
    stimevent.dur(:)=value;
    raw(selected).stimulus(stimnames{idx})=stimevent;
    
    set(handles.figure1,'UserData',raw);
   updateDraw(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function uimenu_setamp_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_setamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');
stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');
stimevent = raw(selected).stimulus(stimnames{idx});

ANSWER = inputdlg('Enter amplutude','Amplitude',1,{'1'});
if(~isempty(ANSWER) && ~isempty(str2num(ANSWER{1})))
    value=str2num(ANSWER{1});
    
    
    stimevent.amp(:)=value;
    raw(selected).stimulus(stimnames{idx})=stimevent;
    
    set(handles.figure1,'UserData',raw);
    updateDraw(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function uimenushift_Callback(hObject, eventdata, handles)
% hObject    handle to uimenushift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');
stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');
stimevent = raw(selected).stimulus(stimnames{idx});

ANSWER = inputdlg('Enter time shift','shift',1,{'0'});
if(~isempty(ANSWER) && ~isempty(str2num(ANSWER{1})))
    value=str2num(ANSWER{1});
    
    
    stimevent.onset=stimevent.onset+value;
    raw(selected).stimulus(stimnames{idx})=stimevent;
    
    set(handles.figure1,'UserData',raw);
    updateDraw(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function uimenu_rename_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');
stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');

job=nirs.modules.RenameStims;
job.listOfChanges{1,1}=stimnames{idx};

prompt={'Enter new stimulus name:'};
name='Rename';
numlines=1;
defaultanswer={stimnames{idx}};

answer=inputdlg(prompt,name,numlines,defaultanswer);
if(~isempty(answer))
    job.listOfChanges{1,2}=answer{1};
    raw(selected)=job.run(raw(selected));
    
    if(idx>length(nirs.getStimNames(raw(selected))))
        set(handles.listbox1,'Value',1);
    end
    
    set(handles.figure1,'UserData',raw);
    updateDraw(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function uimenu_remove_Callback(hObject, eventdata, handles)
% hObject    handle to uimenu_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');
stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');


ButtonName = questdlg(['Remove ' stimnames{idx} '?'], 'Confirm', 'No','Yes','No');
if(strcmp(ButtonName,'Yes'))
    job=nirs.modules.DiscardStims;
    job.listOfStims={stimnames{idx}};
    raw(selected)=job.run(raw(selected));
    
    if(idx>length(nirs.getStimNames(raw(selected))))
        set(handles.listbox1,'Value',1);
    end
    
    set(handles.figure1,'UserData',raw);
    updateDraw(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function stimtutil_Callback(hObject, eventdata, handles)
% hObject    handle to stimtutil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function tablemenu_Callback(hObject, eventdata, handles)
% hObject    handle to tablemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');
stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');

sv=nirs.design.StimulusEvents;
C=get(handles.uitable1,'Data');

if(isempty(C))
    C=cell(1,4); 
end

d=eventdata.EditData;
if(~isempty(str2num(d))) d=str2num(d); end;
C{eventdata.Indices(1),eventdata.Indices(2)}=d;
for i=1:size(C,1)
    if(~isempty(C{i,2}) & ...
            ~isempty(C{i,3}) & ~isempty(C{i,4}))
        
        if(isempty(stimnames))
           if(isempty(C{i,1}))
               C{i,1}='Event-1';
           end
           stimnames={C{i,1}};
           idx=1;
           set(handles.listbox1,'String',stimnames);
          set(handles.listbox1,'value',idx);
           
        end
        
        if(isempty(C{i,1}) || strcmp(C{i,1},stimnames{idx}))
            sv.name=stimnames{idx};
            
            if(isstr( C{i,2}))
                C{i,2}=num2str(str2num(C{i,2}));
            else
                 C{i,2}=num2str(C{i,2});
            end
            
            
            sv.onset(end+1)=str2num(C{i,2});
            sv.amp(end+1)=C{i,3};
            
           if(isstr( C{i,4}))
             C{i,4}=num2str(str2num(C{i,4}));
           else
                C{i,4}=num2str(C{i,4});
           end
            
            sv.dur(end+1)=str2num(C{i,4});
        elseif(~isempty(C{i,1}))
            if(ismember(C{i,1},raw(selected).stimulus.keys))
                sv2=raw(selected).stimulus(C{i,1});
            else
                sv2=nirs.design.StimulusEvents;
            end
            sv2.name=C{i,1};
            
            if(isstr( C{i,2}))
                C{i,2}=num2str(str2num(C{i,2}));
            else
                 C{i,2}=num2str(C{i,2});
            end
            
            
            
            sv2.onset(end+1)=str2num(C{i,2});
            sv2.amp(end+1)=C{i,3};
            
            if(isstr( C{i,4}))
                C{i,4}=num2str(str2num(C{i,4}));
            else
                C{i,4}=num2str(C{i,4});
            end
            sv2.dur(end+1)=str2num(C{i,4});
            raw(selected).stimulus(C{i,1})=sv2;
            
        end
  
    end
end
set(handles.uitable1,'Data',C);
raw(selected).stimulus(stimnames{idx})=sv;

set(handles.figure1,'UserData',raw);
updateDraw(hObject, eventdata, handles);
set(handles.uitable1,'Data',C);
return


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.popupmenu1);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.popupmenu1);
closereq;


% --------------------------------------------------------------------
function uimenuaddrow_Callback(hObject, eventdata, handles)
% hObject    handle to uimenuaddrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=guihandles(gcf);

C=get(handles.uitable1,'Data');
C{end+1,1}=C{1,1};
C{end,2}='0';
C{end,3}=1;
C{end,4}='1';

set(handles.uitable1,'Data',C);


% --------------------------------------------------------------------
function showStimMarks_Callback(hObject, eventdata, handles)
% hObject    handle to showStimMarks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(handles.showStimMarks,'checked'),'on'))
    set(handles.showStimMarks,'checked','off');
else
     set(handles.showStimMarks,'checked','on');
end

updateDraw(hObject, eventdata, handles);
return

% --------------------------------------------------------------------
function showfNIRS_Callback(hObject, eventdata, handles)
% hObject    handle to showfNIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(handles.showfNIRS,'checked'),'on'))
    set(handles.showfNIRS,'checked','off');
else
     set(handles.showfNIRS,'checked','on');
end
updateDraw(hObject, eventdata, handles);
return

% --------------------------------------------------------------------
function showAuxData_Callback(hObject, eventdata, handles)
% hObject    handle to showAuxData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(handles.showAuxData,'checked'),'on'))
    set(handles.showAuxData,'checked','off');
else
     set(handles.showAuxData,'checked','on');
end

updateDraw(hObject, eventdata, handles);
return

% --------------------------------------------------------------------
function click_window_Callback(hObject, eventdata, handles)
% hObject    handle to click_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


h = imrect;
position = wait(h);

onset=position(1);
dur = position(3);

C=get(handles.uitable1,'Data');
if(~isempty(C))
    C{end+1,1}=C{1,1};
else
    C{end+1,1}='Event-1';
end
C{end,2}=num2str(onset);
C{end,3}=1;
C{end,4}=num2str(dur);
set(handles.uitable1,'Data',C);

delete(h);
uitable1_CellEditCallback(hObject, eventdata, handles);

return


% --------------------------------------------------------------------
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function import_Callback(hObject, eventdata, handles)
% hObject    handle to import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function importcurrent_xls_Callback(hObject, eventdata, handles)
% hObject    handle to importcurrent_xls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');

[filename, pathname] = uigetfile( ...
    {'*.xls'}, ...
    'Select file');

if isequal(filename,0) || isequal(pathname,0)
    return
end

filename=fullfile(pathname,filename);
raw(selected)=nirs.design.read_excel2stim(raw(selected),filename);
set(handles.figure1,'UserData',raw);
 StimUtil_GUI('updateDraw');


% --------------------------------------------------------------------
function importall_xls_Callback(hObject, eventdata, handles)
% hObject    handle to importall_xls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


raw=get(handles.figure1,'UserData');

[filename, pathname] = uigetfile( ...
    {'*.xls'}, ...
    'Select file');

if isequal(filename,0) || isequal(pathname,0)
    return
end

filename=fullfile(pathname,filename);
raw=nirs.design.read_excel2stim(raw,filename);
set(handles.figure1,'UserData',raw);
StimUtil_GUI('updateDraw');



% --------------------------------------------------------------------
function import_clipboard_Callback(hObject, eventdata, handles)
% hObject    handle to import_clipboard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data=clipboard('paste');
header=textscan(data(1:min(find(double(data)==13))),'%s','Delimiter','\t');
C={}; cnt=1;
s=[];
for i=1:length(header{1})
    s=[s '%s'];
end
a=textscan(data,s,'Delimiter','\t');

for i=1:length(header{1})-1
    names{i,1}=header{1}{i}(1:max(strfind(header{1}{i},'_'))-1);
end
unames=unique(names);

C={};

for i=1:length(unames)
    lstOn=find(ismember(header{1},[unames{i} '_onset']));
    lstDur=find(ismember(header{1},[unames{i} '_duration']));
    lstAmp=find(ismember(header{1},[unames{i} '_amplitude']));
   
    for j=3:length(a{lstOn})
        C{end+1,1}=unames{i};
        C{end,2}=a{lstOn}{j};
        C{end,3}=str2num(a{lstAmp}{j});
        C{end,4}=a{lstDur}{j};
    end
    
end
    
set(handles.uitable1,'Data',C);
uitable1_CellEditCallback(hObject, eventdata, handles);


return


% --------------------------------------------------------------------
function exportcurrent_xls_Callback(hObject, eventdata, handles)
% hObject    handle to exportcurrent_xls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');

[filename, pathname] = uiputfile( ...
    {'*.xls'}, ...
    'Save as');

if isequal(filename,0) || isequal(pathname,0)
    return
end

filename=fullfile(pathname,filename);
nirs.design.save_stim2excel(raw(selected),filename);
      
    
return


% --------------------------------------------------------------------
function exportall_xls_Callback(hObject, eventdata, handles)
% hObject    handle to exportall_xls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

raw=get(handles.figure1,'UserData');

[filename, pathname] = uiputfile( ...
    {'*.xls'}, ...
    'Save as');

if isequal(filename,0) || isequal(pathname,0)
    return
end

filename=fullfile(pathname,filename);
nirs.design.save_stim2excel(raw,filename);


% --------------------------------------------------------------------
function export_clipboard_Callback(hObject, eventdata, handles)
% hObject    handle to export_clipboard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');


keys=raw(selected).stimulus.keys;
Array={};
for j=1:length(keys)
    
    st=raw(selected).stimulus(keys{j});
    if(isa(st,'nirs.design.StimulusEvents'))
        Array{1,j+3*(j-1)}=[keys{j} '_onset'];
        Array{1,j+1+3*(j-1)}=[keys{j} '_duration'];
        Array{1,j+2+3*(j-1)}=[keys{j} '_amplitude'];
        
        for k=1:length(st.onset)
            Array{1+k,j+3*(j-1)}=st.onset(k);
            Array{1+k,j+1+3*(j-1)}=st.dur(k);
            Array{1+k,j+2+3*(j-1)}=st.amp(k);
        end
    elseif(isa(st,'nirs.design.StimulusVector'))
        Array{1,j+3*(j-1)}=[keys{j} '_vector'];
        Array{1,j+1+3*(j-1)}=[keys{j} '_time'];
        
        for k=1:length(st.vector)
            Array{1+k,j+3*(j-1)}=st.vector(k);
            Array{1+k,j+1+3*(j-1)}=st.time(k);
        end
        
    end
    
    
end
for j=1:size(Array,2)
    
    for k=1:size(Array,1)
        if(isempty(Array{k,j}))
            Array{k,j}=NaN;
        end
    end
end

nirs.util.copytable2clip(cell2table(Array(2:end,:),'VariableNames',Array(1,:)));
disp('data copied to clipboard');

return

% --------------------------------------------------------------------
function ViewData_Callback(hObject, eventdata, handles)
% hObject    handle to ViewData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_removegaps_Callback(hObject, eventdata, handles)
% hObject    handle to menu_removegaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 prompt={'Max Gap duration','Min inter-stim inteval'};
   name='Remove Gaps';
   numlines=1;
   defaultanswer={'2','NaN'};
   
   answer=inputdlg(prompt,name,numlines,defaultanswer);
   job=nirs.modules.RemoveStimGaps;
   job.maxDuration=str2num(answer{1});
   job.minISI=str2num(answer{2});
   job.seperate_conditions=true;
   
   raw=get(handles.figure1,'UserData');
   selected=get(handles.popupmenu1,'value');
   
   stimnames=nirs.getStimNames(raw(selected));
   idx=get(handles.listbox1,'Value');
   
   st=raw(selected).stimulus;
   raw(selected).stimulus=Dictionary;
   raw(selected).stimulus(stimnames{idx})=st(stimnames{idx});
   raw(selected)=job.run(raw(selected));
   st(stimnames{idx})=raw(selected).stimulus(stimnames{idx});
   raw(selected).stimulus=st;
   set(handles.figure1,'UserData',raw);
   StimUtil_GUI('updateDraw');

return

% --------------------------------------------------------------------
function menu_importtemplate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_importtemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_duplicate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_duplicate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');

stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');

raw(selected).stimulus([stimnames{idx} '_copy'])=raw(selected).stimulus(stimnames{idx});
   set(handles.figure1,'UserData',raw);
   StimUtil_GUI('updateDraw');

return

% --------------------------------------------------------------------
function menu_import_Callback(hObject, eventdata, handles)
% hObject    handle to menu_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
return

% --------------------------------------------------------------------
function menu_removeodds_Callback(hObject, eventdata, handles)
% hObject    handle to menu_removeodds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');

stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');
st=raw(selected).stimulus(stimnames{idx});
st.onset(1:2:end)=[];
st.dur(1:2:end)=[];
st.amp(1:2:end)=[];
raw(selected).stimulus(stimnames{idx})=st;
   set(handles.figure1,'UserData',raw);
   StimUtil_GUI('updateDraw');


return

% --------------------------------------------------------------------
function menu_removeevens_Callback(hObject, eventdata, handles)
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');

stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');
st=raw(selected).stimulus(stimnames{idx});
st.onset(2:2:end)=[];
st.dur(2:2:end)=[];
st.amp(2:2:end)=[];
raw(selected).stimulus(stimnames{idx})=st;
   set(handles.figure1,'UserData',raw);
   StimUtil_GUI('updateDraw');
return

% --------------------------------------------------------------------
function menu_removecustom_Callback(hObject, eventdata, handles)
% hObject    handle to menu_removecustom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');


 prompt={'Enter list to remove'};
   name='Remove events';
   numlines=1;
   defaultanswer={'[1:3:end]'};
answer=inputdlg(prompt,name,numlines,defaultanswer);


stimnames=nirs.getStimNames(raw(selected));
idx=get(handles.listbox1,'Value');
st=raw(selected).stimulus(stimnames{idx});
eval(['st.onset(' answer{1} ')=[];']);
eval(['st.dur(' answer{1} ')=[];']);
eval(['st.amp(' answer{1} ')=[];']);

raw(selected).stimulus(stimnames{idx})=st;
   set(handles.figure1,'UserData',raw);
   StimUtil_GUI('updateDraw');

return


% --------------------------------------------------------------------
function addnew_Callback(hObject, eventdata, handles)
% hObject    handle to addnew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');
st=nirs.design.StimulusEvents;
raw(selected).stimulus('New Event')=st;
set(handles.figure1,'UserData',raw);
StimUtil_GUI('updateDraw');

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

% Last Modified by GUIDE v2.5 20-Sep-2017 21:02:42

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

Filenames={};
for i=1:length(raw)
    [~,Filenames{i}]=fileparts(raw(i).description);
end
set(handles.popupmenu1,'String',Filenames);
set(handles.popupmenu1,'Value',1);
set(handles.figure1,'UserData',raw);

StimUtil_GUI('updateDraw');
uiwait(handles.figure1);
return

function updateDraw(varargin)

handles=guihandles(gcf);
raw=get(handles.figure1,'UserData');
selected=get(handles.popupmenu1,'value');

stimnames=nirs.getStimNames(raw(selected));
if(get(handles.listbox1,'Value')<length(stimnames))
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
    C{i,4}=stimevent.dur(i);
end
set(handles.uitable1,'Data',C);


axes(handles.axes1);
cla(handles.axes1);
raw(selected).draw([]);





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
closereq;

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
StimUtil_GUI('updateDraw');


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
StimUtil_GUI('updateDraw');


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
    StimUtil_GUI('updateDraw');
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
    StimUtil_GUI('updateDraw');
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
    StimUtil_GUI('updateDraw');
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
    StimUtil_GUI('updateDraw');
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
    StimUtil_GUI('updateDraw');
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

for i=1:size(C,1)
    if(~isempty(C{i,1} & ~isempty(C{i,2}) & ...
            ~isempty(C{i,3}) & ~isempty(C{i,4})))
        
        if(strcmp(C{i,1},stimnames{selected}))
            sv.name=stimnames{selected};
            sv.onset(end+1)=C{i,2};
            sv.amp(end+1)=C{i,3};
            sv.dur(end+1)=C{i,4};
        elseif(~isempty(C{i,1}))
            if(ismember(C{i,1},raw(selected).stimulus.keys))
                sv2=raw(selected).stimulus(C{i,1});
            else
                sv2=nirs.design.StimulusEvents;
            end
            sv2.name=C{i,1};
            sv2.onset(end+1)=C{i,2};
            sv2.amp(end+1)=C{i,3};
            sv2.dur(end+1)=C{i,4};
            raw(selected).stimulus(C{i,1})=sv2;
            
        end
  
    end
end
raw(selected).stimulus(stimnames{idx})=sv;

set(handles.figure1,'UserData',raw);
StimUtil_GUI('updateDraw');


return


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;
%closereq;

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closereq;


% --------------------------------------------------------------------
function uimenuaddrow_Callback(hObject, eventdata, handles)
% hObject    handle to uimenuaddrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=guihandles(gcf);

C=get(handles.uitable1,'Data');
C{end+1,1}=C{1,1};
C{end,2}=0;
C{end,3}=1;
C{end,4}=1;

set(handles.uitable1,'Data',C);

function varargout = nirsIO_GUI(varargin)
% NIRSIO_GUI MATLAB code for nirsIO_GUI.fig
%      NIRSIO_GUI, by itself, creates a new NIRSIO_GUI or raises the existing
%      singleton*.
%
%      H = NIRSIO_GUI returns the handle to a new NIRSIO_GUI or the handle to
%      the existing singleton*.
%
%      NIRSIO_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NIRSIO_GUI.M with the given input arguments.
%
%      NIRSIO_GUI('Property','Value',...) creates a new NIRSIO_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nirsIO_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nirsIO_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nirsIO_GUI

% Last Modified by GUIDE v2.5 15-Aug-2015 21:35:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nirsIO_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @nirsIO_GUI_OutputFcn, ...
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


% --- Executes just before nirsIO_GUI is made visible.
function nirsIO_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nirsIO_GUI (see VARARGIN)

% Choose default command line output for nirsIO_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nirsIO_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

[mtree, container] = uitree('v0', 'Root',[pwd filesep],...
    'ExpandFcn', @myExpfcn,'Parent',handles.uipanel1); % Parent is ignored
set(container, 'Parent', handles.uipanel1);
mtree.setMultipleSelectionEnabled(true);
set(container,'units','normalized','position',[0 0 1 1]);
set(container,'Tag','IOtree')
set(container,'Userdata',mtree);
set(handles.edit1,'String',[pwd filesep]);
return


% ---------------------------------------------
function nodes = myExpfcn(tree, value)

try
    count = 0;
    ch = dir(value);
    
    for i=1:length(ch)
        if (any(strcmp(ch(i).name, {'.', '..', ''})) == 0)
            
            if ch(i).isdir
                count = count + 1;
                iconpath = [matlabroot, '/toolbox/matlab/icons/foldericon.gif'];
                nodes(count) = uitreenode('v0',[value, ch(i).name, filesep], ...
                    ch(i).name, iconpath, ~ch(i).isdir);
            elseif(~isempty(strfind(ch(i).name,'.nirs')))
                count = count + 1;
                iconpath = [matlabroot, '/toolbox/matlab/icons/pageicon.gif'];
                nodes(count) = uitreenode('v0',[value, ch(i).name, filesep], ...
                    ch(i).name, iconpath, ~ch(i).isdir);
            end
            
        end
    end
catch
    error('MyApplication:UnrecognizedNode', ...
        ['The uitree node type is not recognized. You may need to ', ...
        'define an ExpandFcn for the nodes.']);
end

if (count == 0)
    nodes = [];
end
return
% ---------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = nirsIO_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_changeroot.
function pushbutton_changeroot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_changeroot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
directoryname = uigetdir(pwd, 'Pick a Directory');

directoryname=[directoryname filesep];
delete(findobj('tag','IOtree'))

[mtree, container] = uitree('v0', 'Root',directoryname,...
    'ExpandFcn', @myExpfcn,'Parent',handles.uipanel1); % Parent is ignored
set(container, 'Parent', handles.uipanel1);
mtree.setMultipleSelectionEnabled(true);
set(container,'units','normalized','position',[0 0 1 1]);
set(container,'Tag','IOtree')
set(container,'Userdata',mtree);
set(handles.edit1,'String',directoryname);

return

% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileext='.nirs';

tempdata=get(handles.figure1,'Userdata');
if(isempty(tempdata)); tempdata=nirs.core.Data.empty; end;

mtree=get(findobj('tag','IOtree'),'Userdata');
nodes=mtree.getSelectedNodes;

seldir={}; selfile={};
for idx=1:length(nodes);
    sel=nodes(idx).getValue;
    if(sel(end)==filesep), sel(end)=[]; end;
    if(isdir(sel))
        seldir{end+1}=[sel filesep];
    else
        selfile{end+1}=sel;
    end
end

if(~isempty(selfile))
    d=nirs.io.loadDotNirs(selfile);
    for idx=1:length(d);
        rsplit = strsplit( selfile{idx}, filesep );
        d(idx).demographics('group')=rsplit{end-2};
        d(idx).demographics('subject')= rsplit{end-1};
    end
    tempdata=[tempdata; d];
end

%if selected folder
for idx=1:length(seldir)
    if(seldir{idx}(end)==filesep), seldir{idx}(end)=[]; end;
    files = rdir([seldir{idx} filesep '**' filesep '*' fileext]);
    
    for iFile = 1:length( files )
        fsplit = strsplit( files(iFile).name, filesep );
        rsplit = strsplit( seldir{idx}, filesep );
        demo(iFile) = length(fsplit(length(rsplit)+1:end-1));
    end
    if(any(demo~=mean(demo)));
        warndlg('Files need to have the same folder heiraricah structure');
        return
    end
    demo=demo(1);
    
    switch(demo)
        case 0
            demolabels={};
        case 1
            demolabels={'subject'};
        case 2
            demolabels={'group','subject'};
        case 3
            demolabels={'group','subject','session'};
        otherwise
            demolabels={'group','subject'};
    end
    disp(['Loading a total of ' num2str(length(files)) ' files']);
    d=nirs.io.loadDirectory(seldir{idx},demolabels);
    tempdata=[tempdata; d];
end

set(handles.figure1,'Userdata',tempdata);
updateIOtable;

return

function updateIOtable
handles=guihandles(findobj('tag','figure1','name','nirsIO_GUI'));
tempdata=get(handles.figure1,'Userdata');
demo=nirs.createDemographicsTable(tempdata);
demo=[table({tempdata.description}','VariableNames',{'File_Name'}) demo];

set(handles.uitable1,'Data',table2cell(demo));
set(handles.uitable1,'ColumnName',demo.Properties.VariableNames);

return


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

tempdata=get(handles.figure1,'Userdata');
fileIdx=eventdata.Indices(1);
VarIdx=eventdata.Indices(2);

if(VarIdx==1)
    tempdata(fileIdx).description=eventdata.NewData;
else
    keys=tempdata(fileIdx).demographics.keys;
    tempdata(fileIdx).demographics(keys{VarIdx-1})=eventdata.NewData;
end
set(handles.figure1,'Userdata',tempdata);
updateIOtable;


return


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
% tempdata=get(handles.figure1,'Userdata');
% fileIdx=eventdata.Indices(1);
% disp(tempdata(fileIdx));
return


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tempdata=get(handles.figure1,'Userdata');
name=get(handles.edit_namedata,'string');

d=[];
try; d=evalin('base',name); end

if(~isempty(d))
    ButtonName = questdlg(['Variable ' name ' already exisits in workspace'], ...
        'Warning: Overwrite Data', ...
        'Overwrite', 'Concatinate', 'Cancel', 'Cancel');
    if(strcmp( ButtonName,'Cancel')), return; end;
    if(strcmp( ButtonName,'Concatinate')), tempdata=[d; tempdata]; end;
end
disp(['Creating variable ' name ' in base workspace'])
assignin('base',name,tempdata);


function edit_namedata_Callback(hObject, eventdata, handles)
% hObject    handle to edit_namedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_namedata as text
%        str2double(get(hObject,'String')) returns contents of edit_namedata as a double


% --- Executes during object creation, after setting all properties.
function edit_namedata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_namedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = QCviewer(varargin)
% QCVIEWER MATLAB code for QCviewer.fig
%      QCVIEWER, by itself, creates a new QCVIEWER or raises the existing
%      singleton*.
%
%      H = QCVIEWER returns the handle to a new QCVIEWER or the handle to
%      the existing singleton*.
%
%      QCVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QCVIEWER.M with the given input arguments.
%
%      QCVIEWER('Property','Value',...) creates a new QCVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QCviewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QCviewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QCviewer

% Last Modified by GUIDE v2.5 04-Feb-2020 13:19:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QCviewer_OpeningFcn, ...
                   'gui_OutputFcn',  @QCviewer_OutputFcn, ...
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


% --- Executes just before QCviewer is made visible.
function QCviewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QCviewer (see VARARGIN)

% Choose default command line output for QCviewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QCviewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.figure1,'Userdata',varargin{1});

updateDrawing(handles);

return


function updateDrawing(handles)

chanLst=[1 2 3];  % for now

data = get(handles.figure1,'Userdata');
axes(handles.axes1);
cla(handles.axes1);
h=data.draw(chanLst,false,handles.axes1);
legend off;

lb=get(handles.slider2,'Value');
ub=get(handles.slider3,'Value');

u=data.time(1)+max([lb ub])/100*(data.time(end)-data.time(1));
l=data.time(1)+min([lb ub])/100*(data.time(end)-data.time(1));

a=min(get(handles.axes1,'YLim'));
b=max(get(handles.axes1,'YLim'));
hold on;
c=fill([u l l u],[a a b b],'r');
set(c,'facealpha',.1);
line([u u],[a b],'color','g');
line([l l],[a b],'color','g');

lst=find(data.time>=l & data.time<=u);
cla(handles.axes3); 
axes(handles.axes3);
hold on;
for i=1:length(chanLst)
    [p,f]=pwelch(data.data(lst,chanLst(i))-mean(data.data(lst,chanLst(i))),[],[],[],data.Fs);
    
    plot(handles.axes3,f,log10(p),'color',h(i).Color);
end

d=data.data(lst,chanLst);

tbl=table;
tbl.source=data.probe.link(chanLst,:).source;
tbl.detector=data.probe.link(chanLst,:).detector;
tbl.type=data.probe.link(chanLst,:).type;
tbl.SNI=nirs.math.structnoiseindex(d)';
tbl.Mean=mean(d,1)';
tbl.StdDev=std(d)';
tbl.Skewness=skewness(d)';
tbl.Kurtosis=kurtosis(d)';

for i=1:length(chanLst)
    [~,~,LMC(i,1)]=lmctest(d(:,i));
    [~,~,KPSS(i,1)]=kpsstest(d(:,i));
    [~,~,ROOT(i,1)]=adftest(d(:,i));
end
tbl.LMC_Stationarity=LMC;
tbl.KPSS=KPSS;
tbl.UnitRoot_Stationarity=ROOT;


set(handles.uitable2,'Data',table2cell(tbl));
set(handles.uitable2,'ColumnName',tbl.Properties.VariableNames);


return


% --- Outputs from this function are returned to the command line.
function varargout = QCviewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


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


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

updateDrawing(handles);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

updateDrawing(handles);



% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

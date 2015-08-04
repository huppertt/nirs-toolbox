function varargout = closeSave(varargin)
% CLOSESAVE M-file for closeSave.fig
%      CLOSESAVE, by itself, creates a new CLOSESAVE or raises the existing
%      singleton*.
%
%      H = CLOSESAVE returns the handle to a new CLOSESAVE or the handle to
%      the existing singleton*.
%
%      CLOSESAVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLOSESAVE.M with the given input arguments.
%
%      CLOSESAVE('Property','Value',...) creates a new CLOSESAVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before closeSave_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to closeSave_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help closeSave

% Last Modified by GUIDE v2.5 18-Jun-2010 14:58:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @closeSave_OpeningFcn, ...
                   'gui_OutputFcn',  @closeSave_OutputFcn, ...
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


% --- Executes just before closeSave is made visible.
function closeSave_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to closeSave (see VARARGIN)

% Choose default command line output for closeSave
handles.output = hObject;

handles.SW = varargin{1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes closeSave wait for user response (see UIRESUME)
% uiwait(handles.closeSave);


% --- Outputs from this function are returned to the command line.
function varargout = closeSave_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbYes.
function pbYes_Callback(hObject, eventdata, handles)
% hObject    handle to pbYes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global cSave
cSave = true;

delete(handles.closeSave);

% --- Executes on button press in pbNo.
function pbNo_Callback(hObject, eventdata, handles)
% hObject    handle to pbNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

closeSave_CloseRequestFcn(hObject, eventdata, handles)

% --- Executes when user attempts to close closeSave.
function closeSave_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to closeSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global cSave
cSave = false;

delete(handles.closeSave);

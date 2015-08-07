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

% Last Modified by GUIDE v2.5 06-Aug-2015 21:10:48

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


% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
    tempdata=[tempdata d];
end

%if selected folder
for idx=1:length(seldir)
    if(~isempty(dir(fullfile(seldir{idx},'*.nirs'))))
        % This is a subjects level folder
        d=nirs.io.loadDirectory(seldir{idx},{});
        rsplit = strsplit( seldir{idx}, filesep );
        for idx=1:length(d);
            d(idx).demographics('group')=rsplit{end-2};
            d(idx).demographics('subject')= rsplit{end-1};
        end
    else
        d=nirs.io.loadDirectory(seldir{idx},{'group','subject'});
    end
end




return
function varargout = StimulusWizardDM(varargin)
% STIMULUSWIZARDDM M-file for StimulusWizardDM.fig
%      STIMULUSWIZARDDM, by itself, creates a new STIMULUSWIZARDDM or raises the existing
%      singleton*.
%
%      H = STIMULUSWIZARDDM returns the handle to a new STIMULUSWIZARDDM or the handle to
%      the existing singleton*.
%
%      STIMULUSWIZARDDM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STIMULUSWIZARDDM.M with the given input arguments.
%
%      STIMULUSWIZARDDM('Property','Value',...) creates a new STIMULUSWIZARDDM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StimulusWizardDM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StimulusWizardDM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StimulusWizardDM

% Last Modified by GUIDE v2.5 18-Jun-2010 16:17:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StimulusWizardDM_OpeningFcn, ...
                   'gui_OutputFcn',  @StimulusWizardDM_OutputFcn, ...
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


% --- Executes just before StimulusWizardDM is made visible.
function StimulusWizardDM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StimulusWizardDM (see VARARGIN)

% Choose default command line output for StimulusWizardDM
handles.output = hObject;
handles.isData = false;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes StimulusWizardDM wait for user response (see UIRESUME)
% uiwait(handles.StimulusWizardDM);

if(nargin>3)
    menuLoad_Callback(handles.menuLoad, [], handles,varargin{1})

end

return

% --- Outputs from this function are returned to the command line.
function varargout = StimulusWizardDM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuLoad_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to menuLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Determine Working Directory
if(nargin==4)
loc=varargin{1};
else
loc = uigetdir(pwd, 'Load .NIRS from Directory');
end

if ~isequal(loc,0)

    % Clear Previous Data
    cleartable(handles);

    % Determine Files
    files = dir([loc, '/*.nirs']);

    % Read in Files from Working Directory
    nirs(length(files), 1).name = 0;
    for i=1:length(files)
        
        % Initialize
        clear Include;
        
        % Load Data
        nirs(i).name = files(i).name;
        load([loc, '/', nirs(i).name], '', '-MAT');
        
        if(exist('aux10','var'))
            aux=aux10;
        end
        if(~exist('aux','var'))
            aux=[];
        end
        if(~exist('stiminfo','var'))
            stiminfo=[];
        end
        
        aux(:,end+[1:size(s,2)])=s*750;
        % Save Default Data
        nirs(i).SD = SD;
        nirs(i).aux = aux;
        nirs(i).d = d;
        nirs(i).ml = ml;
        nirs(i).s = s;
        nirs(i).stiminfo = stiminfo;
        nirs(i).t = t;
        
        % Check if File has been Run Before
        if (exist('Include'))
            try
            nirs(i).Include = Include;
            nirs(i).StimDesign = StimDesign;
            nirs(i).StimTitle = StimTitle;
            nirs(i).Condition = Condition;
            nirs(i).AuxiliaryFile = AuxiliaryFile;
            nirs(i).AuxiliaryPath = AuxiliaryPath;
            nirs(i).FastScanFile = FastScanFile;
            nirs(i).FastScanPath = FastScanPath;
            catch
                 nirs(i).Include = false;
            nirs(i).StimDesign = struct([]);
            nirs(i).StimTitle = '';
            nirs(i).Condition = '';
            nirs(i).AuxiliaryFile = '';
            nirs(i).AuxiliaryPath = '';
            nirs(i).FastScanFile = '';
            nirs(i).FastScanPath = '';
            end
        else
            nirs(i).Include = false;
            nirs(i).StimDesign = struct([]);
            nirs(i).StimTitle = '';
            nirs(i).Condition = '';
            nirs(i).AuxiliaryFile = '';
            nirs(i).AuxiliaryPath = '';
            nirs(i).FastScanFile = '';
            nirs(i).FastScanPath = '';
        end
        
        % Check Comments
        if (exist('Comments'))
            nirs(i).Comments = Comments;
        else
            nirs(i).Comments = '';
        end
        
        % Populate Table
        tbdata{i,1} = nirs(i).name;
        tbdata{i,2} = nirs(i).Include;
        tbdata{i,3} = nirs(i).StimTitle;
        tbdata{i,4} = nirs(i).Condition;
        tbdata{i,5} = nirs(i).AuxiliaryFile;
        tbdata{i,6} = nirs(i).FastScanFile;
        tbdata{i,7} = nirs(i).Comments;
        
    end
 
    % Save Data
    handles.loc = loc;
    handles.files = files;
    handles.nirs = nirs;
    handles.tbdata = tbdata;
    handles.isData = true;

    guidata(hObject, handles);
    updatetable(handles);

end
    

% --------------------------------------------------------------------
function menuSave_Callback(hObject, eventdata, handles)
% hObject    handle to menuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable,'enable','on');

if (handles.isData)
    for i=1:length(handles.nirs)

        if (handles.nirs(i).Include)
            tempstruct = handles.nirs(i);
        else
            tempstruct = handles.nirs(i);
            tempstruct.StimDesign = struct([]);
            tempstruct.StimTitle = '';
            tempstruct.Condition = '';
            tempstruct.AuxiliaryFile = '';
            tempstruct.AuxiliaryPath = '';
            tempstruct.FastScanFile = '';
            tempstruct.FastScanPath = '';
            % tempstruct.Comments = ''; % -- Always Save Comments
        end
        save([handles.loc, '/', handles.nirs(i).name], '-STRUCT', 'tempstruct');

    end
    
    msgbox('NIRS Saved');
end

% --------------------------------------------------------------------
function menuSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to menuSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.isData)
    
    [lead, path] = uiputfile({'*.nirs', 'NIRS Files (*.nirs)'}, 'Save As (File Header)', 'SW_');
    p = find(lead=='.');
    if ~isempty(p)
        lead = lead(1:p(1)-1);
    end
    
    for i=1:length(handles.nirs)
        
        if (handles.nirs(i).Include)
            tempstruct = handles.nirs(i);
        else
            tempstruct = handles.nirs(i);
            tempstruct.StimDesign = struct([]);
            tempstruct.StimTitle = '';
            tempstruct.Condition = '';
            tempstruct.AuxiliaryFile = '';
            tempstruct.AuxiliaryPath = '';
            tempstruct.FastScanFile = '';
            tempstruct.FastScanPath = '';
            tempstruct.Comments = '';
        end
        save([path, '/', [lead, handles.nirs(i).name]], '-STRUCT', 'tempstruct');
        
    end
    
    msgbox('NIRS Saved');
end

% --------------------------------------------------------------------
function menuExit_Callback(hObject, eventdata, handles)
% hObject    handle to menuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
StimulusWizardDM_CloseRequestFcn(handles.StimulusWizardDM, eventdata, handles);


% --- Executes when selected cell(s) is changed in uitable.
function uitable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

if (handles.isData) && (numel(eventdata.Indices)==2)
    r = eventdata.Indices(1);
    c = eventdata.Indices(2);

    if (c == 2) % Include? -------------------------------------
        handles.tbdata{r,c} = ~handles.tbdata{r,c};
        handles.nirs(r).Include = ~handles.nirs(r).Include;
        
    elseif (c == 3) % Select Stimulus File ---------------------
        global SD;
        global ST;
        global aux;
        global taux;
        
        SD = handles.nirs(r).StimDesign;
        ST = handles.nirs(r).StimTitle;
        aux = handles.nirs(r).aux;
        taux = handles.nirs(r).t;
        
        h = getStimDesign;
        set(handles.uitable, 'Enable', 'off');
        
        uiwait(h);
        
        set(handles.uitable, 'Enable', 'on');
        handles.nirs(r).StimDesign = SD;
        handles.nirs(r).StimTitle = ST;
        handles.tbdata{r,c} = ST;
        
    elseif (c == 4) % Write Condition --------------------------
        cond = inputdlg('Condition');
        if ~isempty(cond)
            handles.nirs(r).Condition = cond{1};
            handles.tbdata{r,c} = cond{1};
        end
        
    elseif (c == 5) % Select Auxiliary File --------------------
        [file, path] = uigetfile('*.*', 'Select Auxiliary File');
        
        if ~isequal(file, 0) && ~isequal(path, 0)
            handles.nirs(r).AuxiliaryFile = file;
            handles.nirs(r).AuxiliaryPath = path;
            handles.tbdata{r,c} = file;
        end
    
    elseif (c == 6) % Select Fast Scan File --------------------
        [file, path] = uigetfile('*.*', 'Select Fast Scan File');
        
        if ~isequal(file, 0) && ~isequal(path, 0)
            for i=1:length(handles.nirs)
                handles.nirs(i).FastScanFile = file;
                handles.nirs(i).FastScanPath = path;
                handles.tbdata{i,c} = file;
            end
        end
        
    elseif (c == 7) % Write Comments ---------------------------
        com = inputdlg('Comments');
        if ~isempty(com)
            handles.nirs(r).Comments = com{1};
            handles.tbdata{r,c} = com{1};
        end
        
    end

    guidata(hObject, handles);
    updatetable(handles);
end


% --------------------------------------------------------------------
function cleartable(handles)
% Initialize Table
if handles.isData
    
    rmfield(handles, 'nirs');
    
    for i=1:size(handles.tbdata, 1)
        handles.tbdata{i,1} = '';
        handles.tbdata{i,2} = false;
        handles.tbdata{i,3} = '';
        handles.tbdata{i,4} = '';
        handles.tbdata{i,5} = '';
        handles.tbdata{i,6} = '';
        handles.tbdata{i,7} = '';
    end
    
    updatetable(handles);
    guidata(handles.StimulusWizardDM, handles);
    
end


% --------------------------------------------------------------------
function updatetable(handles)
% Update Table Info
set(handles.uitable, 'Data', handles.tbdata);


% --- Executes when entered data in editable cell(s) in uitable.
function uitable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitable_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close StimulusWizardDM.
function StimulusWizardDM_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to StimulusWizardDM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global cSave;

h = closeSave(hObject);
uiwait(h);

if cSave
    menuSave_Callback(handles.menuSave, eventdata, handles);
end

delete(hObject);

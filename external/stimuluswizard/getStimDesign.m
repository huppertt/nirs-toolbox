function varargout = getStimDesign(varargin)
% GETSTIMDESIGN M-file for getStimDesign.fig
%      GETSTIMDESIGN, by itself, creates a new GETSTIMDESIGN or raises the existing
%      singleton*.
%
%      H = GETSTIMDESIGN returns the handle to a new GETSTIMDESIGN or the handle to
%      the existing singleton*.
%
%      GETSTIMDESIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETSTIMDESIGN.M with the given input arguments.
%
%      GETSTIMDESIGN('Property','Value',...) creates a new GETSTIMDESIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getStimDesign_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to getStimDesign_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getStimDesign

% Last Modified by GUIDE v2.5 20-Dec-2013 11:52:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @getStimDesign_OpeningFcn, ...
                   'gui_OutputFcn',  @getStimDesign_OutputFcn, ...
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

end


% --- Executes just before getStimDesign is made visible.
function getStimDesign_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getStimDesign (see VARARGIN)

% Choose default command line output for getStimDesign
handles.output = hObject;

global SD;
global ST;

% Populate StimDesign Table
tbdata{1,1} = [];
tbdata{1,2} = [];
tbdata{1,3} = [];
tbdata{1,4} = '';
if exist('SD', 'var') && ~isempty(SD)
    tbdata = cell(length(SD), 4);
    for i=1:length(SD)
        tbdata{i,1} = SD(i).onset;
        tbdata{i,2} = SD(i).dur;
        tbdata{i,3} = SD(i).amp;
        tbdata{i,4} = SD(i).cond;
    end 
end
handles.tbdata = tbdata;

% Determine StimDesign Title
contents = get(handles.menuStim,'String');
if exist('ST', 'var') && ~isempty(ST)
    if isequal(ST, contents{1})
        stimtitle = ST;
        set(handles.menuStim, 'Value', 1);
    else
        stimtitle = contents{2};
        set(handles.menuStim, 'Value', 2);
    end
else
    stimtitle = contents{1};
    set(handles.menuStim, 'Value', 1);
end
% if exist('ST', 'var') && ~isempty(ST)
%     for i=1:length(contents)
%         if isequal(contents{i}, ST)
%             set(handles.menuStim, 'Value', i);
%             stimtitle = ST;
%             break;
%         end
%     end
% end
handles.stimtitle = stimtitle;

updatetable(handles);
guidata(hObject, handles);

% UIWAIT makes getStimDesign wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end


% --- Outputs from this function are returned to the command line.
function varargout = getStimDesign_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end


% --- Executes on selection change in menuStim.
function menuStim_Callback(hObject, eventdata, handles)
% hObject    handle to menuStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuStim

ind = get(hObject,'Value');
contents = get(hObject,'String');
handles.stimtitle = contents{ind};

if(ind==3)
[file, path,ind] = uigetfile({'*.txt','*.PDAT'}, 'Select E-Prime File');
if (ind == 2) % E-Prime (SE: *.PDAT)
   [tbdata, run] = readEPrimePDAT(file, path);
   handles.stimtitle = [file, ' Run #', run{1}];
elseif (ind == 1) % E-Prime (MNR: *.txt)
    [tbdata,stimdesign]=readEprimeTXT(fullfile(path,file));
end
end
if(ind==2)
    %Re
    global aux;
    global taux;
    [tbdata]=get_stimdesign_aux(taux,aux);
end
handles.tbdata = tbdata;

updatetable(handles);
guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function menuStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in pbSubmit.
function pbSubmit_Callback(hObject, eventdata, handles)
% hObject    handle to pbSubmit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tbdata = get(handles.uitable, 'Data');
[r, c] = size(handles.tbdata);

global ST;
global SD;

SD = struct([]);

for i=1:r
    SD(i).onset = handles.tbdata{i,1};
    SD(i).dur = handles.tbdata{i,2};
    SD(i).amp = handles.tbdata{i,3};
    SD(i).cond = handles.tbdata{i,4};
end 
ST = handles.stimtitle;

guidata(hObject, handles);

close(ancestor(hObject,'figure'));

end

% --- Executes on button press in pbCancel.
function pbCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pbCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(ancestor(hObject,'figure'));

end

% --- Executes on button press in pbGraph.
function pbGraph_Callback(hObject, eventdata, handles)
% hObject    handle to pbGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.uitable, 'Enable', 'off');

% Enter Auxiliary Channel Data
global aux;
global taux;

% Update Table
handles.tbdata = get(handles.uitable, 'Data');
if ~isfield(handles, 'tbdata_old')
    handles.tbdata_old = handles.tbdata;
end

% Initialize Figure
fig = 10;
if ishandle(fig)
    close(fig);
end
figure(fig);
set(fig, 'Toolbar', 'figure');
pos = get(fig, 'OuterPosition');
set(gcf, 'OuterPosition', [pos(1) pos(2) pos(3)*2 pos(4)*2]);

% Plot AUX vs TIME
naux = size(aux,2) ;
h1 = subplot(2,1,1);
hold all;
leg = cell(1, naux);

for i=1:naux
    plot(taux, aux(:,i));
    leg{i} = ['AUX #', num2str(i)];
end

xlabel('TIME (SEC)'); ylabel('AMPLITUDE');
xlim([min(taux), max(taux)]);
mxAux = max(max(aux(:,1:naux)));
mnAux = min(min(aux(:,1:naux)));
bfAux = ((mxAux + mnAux)/2)*0.1;
ylim([mnAux-bfAux, mxAux+bfAux]);
legend(leg, 'Location', 'SouthOutside');


% Plot Stimulus Conditions
h2 = subplot(2,1,2);
hold all;

cellemp = cellfun(@isempty, handles.tbdata);
if isempty(find(cellemp(:,1:3)==1)) %#ok<EFIND>
    con = handles.tbdata(:,4);
    uniq = unique(con);

    j = 0;
    leg = {};
    mnTime = inf;
    mxTime = 0;
    mnStim = inf;
    mxStim = 0;
    for i=1:length(uniq)

        tf = ismember(con, uniq(i));
        numI = sum(tf);
%         if numI > 1
            j = j + 1;
            leg{j} = uniq{i};

            ons = handles.tbdata(tf,1);
            dur = handles.tbdata(tf,2);
            amp = handles.tbdata(tf,3);

            offset = 0;

            x = [];
            y = [];
            p = 0;
            for k=1:numI
                p = p + 1;
                x(p) = ons{k};
                y(p) = 0 + offset;
                p = p + 1;
                x(p) = ons{k};
                y(p) = amp{k} + offset;
                p = p + 1;
                x(p) = ons{k} + dur{k};
                y(p) = amp{k} + offset;
                p = p + 1;
                x(p) = ons{k} + dur{k};
                y(p) = 0 + offset;

            end

            if (min(x) < mnTime)
                mnTime = min(x);
            end
            if (max(x) > mxTime)
                mxTime = max(x);
            end
            if (min(y) < mnStim)
                mnStim = min(y);
            end
            if (max(y) > mxStim)
                mxStim = max(y);
            end

            plot(x, y);

%         end

    end

    if ~isinf(mnTime) && ~isinf(mnStim)
    %     bfTime = (mnTime + mxTime)/2 * 0.1;
        bfStim = (mnStim + mxStim)/2 * 0.1;
        xlim(xlim(h1));
        ylim([mnStim-bfStim, mxStim+bfStim]);
        legend(leg, 'Location', 'SouthOutside');
    end
else
    msgbox('Incomplete Stim Design');
end
xlabel('TIME (SEC)'); ylabel('AMPLITUDE');

% Initialize Push Buttons & Nested Callbacks
htext = uicontrol(fig, 'Style', 'text', 'String', '', 'Position', [10 10 290 15]);
handles.htext = htext;

uicontrol(fig, 'Style', 'pushbutton', 'String', 'Align', 'Position', [10 35 90 40], 'Callback', @(varargin)pbAlign(fig, hObject, eventdata, handles));
function pbAlign(fig, hObject, eventdata, handles)
    h = datacursormode(fig);
    set(h, 'DisplayStyle', 'window', 'SnapToDataVertex', 'on', 'Enable', 'on');
    
    set(handles.htext, 'String', 'Select AUX Point to Align, then Press Enter.');
    pause;
    s1 = getCursorInfo(h);
    
    set(handles.htext, 'String', 'Select Stimulus Point to Align, then Press Enter.');
	pause
    s2 = getCursorInfo(h);
    
    set(h, 'DisplayStyle', 'datatip', 'Enable', 'off');
    set(handles.htext, 'String', 'Simtulus Design Aligned to AUX.');
    
    handles.tbdata_old = handles.tbdata;
    diff = s2.Position(1) - s1.Position(1);
    for i=1:size(handles.tbdata,1)
        handles.tbdata{i,1} = handles.tbdata{i,1} - diff;
    end
    
    updatetable(handles);
    guidata(handles.output, handles);
    
    close(fig);
    pbGraph_Callback(hObject, eventdata, handles);
end

uicontrol(fig, 'Style', 'pushbutton', 'String', 'Accept', 'Position', [110 35 90 40], 'Callback', @(varargin)pbAccept(fig, hObject, eventdata, handles));
function pbAccept(fig, hObject, eventdata, handles)
    close(fig);
    set(handles.uitable, 'Enable', 'on');
end

uicontrol(fig, 'Style', 'pushbutton', 'String', 'Reject', 'Position', [210 35 90 40], 'Callback', @(varargin)pbReject(fig, hObject, eventdata, handles));
function pbReject(fig, hObject, eventdata, handles)
    handles.tbdata = handles.tbdata_old;
    updatetable(handles);
    guidata(handles.output, handles);
    
    close(fig);
    set(handles.uitable, 'Enable', 'on');
%     pbGraph_Callback(hObject, eventdata, handles);
end

guidata(hObject, handles);

end



% --- Executes on button press in pbAdd.
function pbAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pbAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tbdata = get(handles.uitable, 'Data');
[r, c] = size(handles.tbdata);

handles.tbdata{r+1, 1} = [];
handles.tbdata{r+1, 2} = [];
handles.tbdata{r+1, 3} = [];
handles.tbdata{r+1, 4} = '';

updatetable(handles);
guidata(hObject, handles);

end


% --------------------------------------------------------------------
function updatetable(handles)
% Update Table Info
set(handles.uitable, 'Data', handles.tbdata);

end


% --------------------------------------------------------------------
function cleantable(hObject, handles)
% Remove Invalid Entries
data = get(handles.uitable, 'Data');
del = [];
d = 0;

for i=1:size(data,1)
    if isempty(data{i,1}) || isempty(data{i,2}) || ...
       isempty(data{i,3}) || isempty(data{i,4})
        
        d = d + 1;
        del(d,1) = i;
   
    end
end

data(del,:) = [];
handles.tbdata = data;

updatetable(handles);
guidata(hObject, handles);

end


% --- Executes on key press with focus on uitable and none of its controls.
function uitable_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
eventdata.Key;

end


% --- Executes when selected cell(s) is changed in uitable.
function uitable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
eventdata.Indices;

end


% --- Executes on button press in pbInsert.
function pbInsert_Callback(hObject, eventdata, handles)
% hObject    handle to pbInsert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tbdata = get(handles.uitable, 'Data');
[r, c] = size(handles.tbdata);

inp = inputdlg('Insert Row #:');
row = str2double(inp{1});
if (row > r)
    row = r + 1;
elseif (row < 1)
    row = 1;
end

handles.tbdata = [ handles.tbdata(1:row-1,:) ; ...
                   { [] [] [] '' } ; ...
                   handles.tbdata(row:end,:) ];

updatetable(handles);
guidata(hObject, handles);

end

% --- Executes on button press in pbRemove.
function pbRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pbRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tbdata = get(handles.uitable, 'Data');
[r, c] = size(handles.tbdata);

inp = inputdlg('Remove Row #:');
row = str2double(inp{1});

if (row > 0) && (row <= r)
    
    handles.tbdata = [ handles.tbdata(1:row-1,:) ; ...
                       handles.tbdata(row+1:end,:) ];
    
else
    msgbox('Fail');
end

updatetable(handles);
guidata(hObject, handles);

end


% --- Executes on button press in pbCommand.
function pbCommand_Callback(hObject, eventdata, handles)
% hObject    handle to pbCommand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tbdata = get(handles.uitable, 'Data');
[r, c] = size(handles.tbdata);

inp1 = inputdlg('Enter Onset Command:  E.G. 10:20:300');
clear ons;
if ~isempty(inp1)
    ons = eval(inp1{1});
    
    inp2 = inputdlg('Enter Constant Duration:  E.G. 10');
    clear dur;
    if ~isempty(inp2)
        dur = eval(inp2{1});
    end
end

if exist('ons', 'var') && exist('dur', 'var')
    for i=1:length(ons)
        handles.tbdata{r+i,1} = ons(i);
        handles.tbdata{r+i,2} = dur;
        handles.tbdata{r+i,3} = 1;
        handles.tbdata{r+i,4} = 'A';
    end

end
    
updatetable(handles);
guidata(hObject, handles);

end
    


% --------------------------------------------------------------------
function uitable_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


% --------------------------------------------------------------------
function copytable_Callback(hObject, eventdata, handles)
% hObject    handle to copytable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdata=get(handles.uitable,'Data');
copytable2clip(tdata);

end


% --------------------------------------------------------------------
function pastetable_Callback(hObject, eventdata, handles)
% hObject    handle to pastetable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdata=pasteclip2table;
set(handles.uitable,'Data',tdata);


end


% --------------------------------------------------------------------
function CopyPaste_Callback(hObject, eventdata, handles)
% hObject    handle to CopyPaste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

classdef dontshowassistantdlg < handle
  
  properties (Access = public)
    Figure
    DlgChoice
  end
  
  properties (Access = private)
    
    hndl_msg1
    hndl_msg2
    hndl_msg3
    hndl_fileName
    hndl_filePath
    hndl_info
    
    hndl_checkbox
    hndl_button1
    hndl_button2
    
  end
  
  methods    
    
    %----------------------------------------------------------------------
    function obj = dontshowassistantdlg
      % Constructs the Dialog object
      obj.DlgChoice = '';      
    end    
    
    %----------------------------------------------------------------------
    function render(obj)
      % Renders the popup dialog      
            
      uicFontName = get(0,'DefaultUiControlFontName');
      uicFontSize = get(0,'DefaultUiControlFontSize');
      txtFontName = get(0,'DefaultTextFontName');
      txtFontSize = get(0,'DefaultTextFontSize');      
      
      pf = get(0,'ScreenPixelsPerInch')/96;
      if isunix,
        pf = 1;
      end
      
      % Put up the dialog's figure.
      figw = 475;
      figh = 185;
      obj.Figure = figure('MenuBar', 'None', ...
        'Resize', 'Off', ...
        'Tag', 'designFilterDialog', ...
        'Name', getString(message('signal:designfilt:AssistantHeader')),...
        'IntegerHandle', 'Off', ...
        'HandleVisibility', 'Off', ...
        'NumberTitle', 'Off', ...
        'Position', [450 450 figw figh]*pf, ...
        'Color', [1 1 1], ...
        'Visible', 'Off');
      
      fpos = get(obj.Figure,'Position');
          
      lowerpanel = uipanel(obj.Figure,...
        'Units',get(obj.Figure,'Units'), ...
        'Position', [0, 0, fpos(3)+10, 40]*pf);
      
      % Put up the don't show me again check box.
      obj.hndl_checkbox = uicontrol(lowerpanel, ...
        'Style', 'Checkbox', ...
        'Units',get(obj.Figure,'Units'),...
        'Position', [30 12 200 20], ...
        'FontName',uicFontName,...
        'FontSize',uicFontSize,...
        'String', getString(message('signal:designfilt:DoNotShowThisMessageAgain')), ...
        'HorizontalAlignment', 'Left',...
        'Tag', 'DialogCheckBox', ...
        'Callback',{@checkbox_cb, obj});
      
      obj.hndl_button1 = uicontrol(lowerpanel, ...
        'Style', 'PushButton', ...
        'String', getString(message('signal:designfilt:Yes')), ...
        'FontName',uicFontName,...
        'FontSize',uicFontSize,...
        'Units',get(obj.Figure,'Units'),...
        'Position', [270 8 70 25] , ...
        'HorizontalAlignment', 'Center',...
        'Tag','DialogButton1',...
        'Callback',{@button1_cb, obj});
      
      obj.hndl_button2 = uicontrol(lowerpanel, ...
        'Style', 'PushButton', ...
        'String', getString(message('signal:designfilt:No')), ...
        'FontName',uicFontName,...
        'FontSize',uicFontSize,...
        'Units',get(obj.Figure,'Units'),...
        'Position', [370 8 70 25], ...
        'HorizontalAlignment', 'Center',...
        'Tag','DialogButton2',...
        'Callback',{@button2_cb, obj});
      
      obj.hndl_msg1 = uicontrol(obj.Figure, 'Style', 'Text', ...
        'HorizontalAlignment', 'Left', ...
        'Units',get(obj.Figure,'Units'),...
        'FontName',txtFontName, ...
        'FontSize',txtFontSize + 2,...
        'ForegroundColor','b', ...
        'Position', [15,120,450,55],...
        'BackgroundColor',[1 1 1],...
        'Tag','DialogMsg1');
      
      obj.hndl_msg2 = uicontrol(obj.Figure, 'Style', 'Text', ...
        'HorizontalAlignment', 'Left', ...
        'Units',get(obj.Figure,'Units'),...
        'Position', [15,60,425,60],...
        'FontName',txtFontName, ...
        'FontSize',txtFontSize,...
        'ForegroundColor',[0.3 0.3 0.3], ...
        'BackgroundColor',[1 1 1],...
        'Tag','DialogMsg2');
      
      obj.hndl_msg3 = uicontrol(obj.Figure, 'Style', 'Text', ...
        'HorizontalAlignment', 'Left', ...
        'Units',get(obj.Figure,'Units'),...
        'Position', [15,45,425,25],...
        'FontName',txtFontName, ...
        'FontSize',txtFontSize,...
        'FontWeight','bold',...
        'BackgroundColor',[1 1 1],...
        'Tag','DialogMsg3');
      
       obj.hndl_fileName = uicontrol(obj.Figure, 'Style', 'Text', ...
        'HorizontalAlignment', 'Left', ...
        'Units',get(obj.Figure,'Units'),...
        'Position', [400,100,425,25],...
        'FontName',uicFontName, ...
        'FontSize',uicFontSize+1,...
        'FontWeight','bold',...
        'ForegroundColor',[62 81 113]/255, ...
        'BackgroundColor',[1 1 1],...
        'Tag','fileNameText');
      
      
       obj.hndl_filePath = uicontrol(obj.Figure, 'Style', 'Text', ...
        'HorizontalAlignment', 'Left', ...
        'Units',get(obj.Figure,'Units'),...
        'Position', [400,45,425,25],...
        'FontName',txtFontName, ...
        'FontSize',txtFontSize-2,...
        'ForegroundColor',[0.3 0.3 0.3], ...
        'BackgroundColor',[1 1 1],...
        'Tag','filePathText');
      
       obj.hndl_info = uicontrol(obj.Figure, 'Style', 'Text', ...
        'HorizontalAlignment', 'Left', ...
        'Units',get(obj.Figure,'Units'),...
        'Position', [400,45,425,25],...
        'FontName',txtFontName, ...
        'FontSize',txtFontSize,...
        'ForegroundColor',[0.3 0.3 0.3], ...
        'BackgroundColor',[1 1 1],...
        'Tag','infoText');
     
    end    
    
    %----------------------------------------------------------------------
    function setUicontrolStrings(obj,strs)
      
      if ~iscell(strs)
        strs = {strs};
      end
      
      N = numel(strs);
      
      if N == 1
        set(obj.hndl_checkbox,'String',strs{1});
      elseif N == 2
        set(obj.hndl_checkbox,'String',strs{1});
        set(obj.hndl_button2,'String',strs{2});
      else
        set(obj.hndl_checkbox,'String',strs{1});
        set(obj.hndl_button2,'String',strs{2});        
        set(obj.hndl_button1,'String',strs{3});
      end
      
    end
    
    %----------------------------------------------------------------------
    function setDialogStrings(obj,strs,fileNameMode)
      % Sets the Dialog Strings and adjusts the heigth of the dialog
      % depending on the length of the strings passed in.
      
      if nargin < 3
        fileNameMode = false;
      end
      
      if ~fileNameMode
        [msg1Pos,msg2Pos,msg3Pos,fPos] = getMode1DialogHeight(strs);
        
        set(obj.hndl_msg2,'Visible','on');
        set(obj.hndl_fileName,'Visible','off');
        set(obj.hndl_filePath,'Visible','off');
        set(obj.hndl_info,'Visible','off');
        
        set(obj.hndl_msg1,'String',strs{1},'Position',msg1Pos);
        set(obj.hndl_msg2,'String',strs{2},'Position',msg2Pos);
        set(obj.hndl_msg3,'String',strs{3},'Position',msg3Pos);
        
        set(obj.Figure,'Position',fPos);      
        
      else
        
        [msg1Pos,fileNamePos,filePathPos,infoPos,msg3Pos,fPos,strSplit] = getMode2DialogHeight(strs);
      
        set(obj.hndl_msg2,'Visible','off');
        set(obj.hndl_fileName,'Visible','on');
        set(obj.hndl_filePath,'Visible','on');
        set(obj.hndl_info,'Visible','on');        
      
        % Set Strings and Positions of text boxes
        set(obj.hndl_msg1,'String',strs{1},'Position',msg1Pos);
        set(obj.hndl_fileName,'String',strs{2},'Position',fileNamePos);
        set(obj.hndl_filePath,'String',strSplit,'Position',filePathPos);
        set(obj.hndl_info,'String',strs{4},'Position',infoPos);
        set(obj.hndl_msg3,'String',strs{5},'Position',msg3Pos);
        
        set(obj.Figure,'Position',fPos);
        
      end
      
    end    
        
  end
  
end

%--------------------------------------------------------------------------
function checkbox_cb(~,~,obj)
% obj function runs when the app is closed

obj.DlgChoice = 'check';
delete(obj.Figure)

end
%--------------------------------------------------------------------------
function button1_cb(~,~,obj)
% obj function runs when the app is closed

obj.DlgChoice = 'Yes';
delete(obj.Figure)

end

%--------------------------------------------------------------------------
function button2_cb(~,~,obj)
% obj function runs when the app is closed

obj.DlgChoice = 'No';
delete(obj.Figure)

end

%--------------------------------------------------------------------------
function [msg1Pos,msg2Pos,msg3Pos,fPos] = getMode1DialogHeight(strs)

% Set Initial Position of the Message Boxes
msg1Pos = [15,120,450,55];
msg2Pos = [15,60,450,60];
msg3Pos = [15,45,450,25];

% Initialize variables
numStr = 3;                         % Number of string inputs
wrapLineCount = zeros(1,numStr);    % Count of number of lines wrapped
stSize = zeros(1,numStr);           % Number of lines of each input stings
wrapLen = [55,77,64];               % wrap lengths for different fonts

% Loop Through each String
for i = 1:numStr
  
  st = strs{i};                             % Get each String
  nst = regexp(st, '[\f\n\r]', 'split');    % put lines into cells
  stSize(i) = size(nst,2);                  % Number of lines
  
  % Loop Through each line and check if it wraps around
  for j = 1:stSize(i)
    % If line exceeds max chars for line then increment count 
    if length(nst{j})/wrapLen(i) > 1
      wrapLineCount(i) = wrapLineCount(i) + ceil(length(nst{j})/wrapLen(i)) - 1;
    end   
  end    
end

% Determine the Number of lines required to fit the strings
st1H = (stSize(1) + wrapLineCount(1)) * 25;
st2H = (stSize(2) + wrapLineCount(2)) * 16;
st3H = (stSize(3) + wrapLineCount(3)) * 20;

% Adjust ypos and height for message boxes
msg3Pos(2) = 45;
msg3Pos(4) = st3H;

msg2Pos(2) = msg3Pos(2) + msg3Pos(4) + 10;
msg2Pos(4) = st2H;

msg1Pos(2) = msg2Pos(2) + msg2Pos(4) + 7;
msg1Pos(4) = st1H;

% Adjust height for the figure
h = msg1Pos(2) + msg1Pos(4) + 10;
fPos = [450 450 475 h];

end

%--------------------------------------------------------------------------
function [msg1Pos,fileNamePos,filePathPos,infoPos,msg3Pos,fPos,strSplit] = getMode2DialogHeight(strs)

  msg1Pos = [15,120,450,55];
  msg3Pos = [15,45,450,25];
  filePathPos = [15 90 450 25];
  fileNamePos = [15 110 450 20];  
  infoPos = [15 110 450 20];
  
  hsplit = 10;
  
  fileName = strs{3};  
  [h,strSplit] = getFileNameHeight(fileName);  
  
  % Determine Positions of the text boxes
  infoPos(2) = msg3Pos(2) + msg3Pos(4) + hsplit;
  infoPos(4) = 17;
  
  filePathPos(2) = infoPos(2) + infoPos(4) + hsplit;
  filePathPos(4) = h;
  
  fileNamePos(2) = filePathPos(2) + filePathPos(4);  
  fileNamePos(4) = 17;
  
  msg1Pos(2) = fileNamePos(2) + fileNamePos(4) + hsplit;
  msg1Pos(4) = 25;
    
  figureHeight = msg1Pos(2) + msg1Pos(4) + 10;
  fPos = [450 450 475 figureHeight];

end

%--------------------------------------------------------------------------
function [h,strSplit] = getFileNameHeight(fName)

tempstr = fName;
Count = 1;
strSplit = {};

wordLim = 92;
slashType = filesep;

while length(tempstr)> wordLim
  
  slashPos = strfind(tempstr,slashType);

  p = max(slashPos(slashPos<wordLim));
  strSplit = [strSplit tempstr(1:p)]; %#ok<AGROW>
  
  tempstr = tempstr(p+1:end);
  
  Count = Count+1;
end

strSplit = [strSplit tempstr];

h = Count*17;

end





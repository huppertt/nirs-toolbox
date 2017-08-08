function varargout = designfilt(varargin)
%DESIGNFILT Design a digital filter
%   D = DESIGNFILT(RESP,P1,V1,P2,V2,...) designs a digital filter D with
%   response RESP and specifications P1,V1, ... PN,VN. See <a href="matlab:doc designfilt">documentation</a> to
%   learn about valid syntaxes and the available property options for each
%   filter response, or type DESIGNFILT or DESIGNFILT(RESP) to launch a
%   Filter Design Assistant where you can specify your filter and get MATLAB 
%   code on the command line that executes and designs the filter.
%
%   When you specify an incorrect set of name-value pairs, you are given a
%   choice to launch a Filter Design Assistant that helps specify your
%   filter and generates the correct MATLAB code for you.
%
%   DESIGNFILT(RESP) launches a Filter Design Assistant for the specified
%   filter response RESP.
%
%   RESP can be one of the following:
%
%   'arbmagfir'          Arbitrary magnitude FIR response
%   'bandpassfir'        Bandpass FIR response
%   'bandpassiir'        Bandpass IIR response
%   'bandstopfir'        Bandstop FIR response
%   'bandstopiir'        Bandstop IIR response
%   'differentiatorfir'  Differentiator FIR response
%   'highpassfir'        Highpass FIR response
%   'highpassiir'        Highpass IIR response
%   'hilbertfir'         Hilbert transformer FIR response
%   'lowpassfir'         Lowpass FIR response
%   'lowpassiir'         Lowpass IIR response
%
%   DESIGNFILT with no inputs launches a response selector where you choose
%   the filter response RESP used to design your filter.
%
%   DESIGNFILT(D) launches a Filter Design Assistant so that you can
%   edit the digital filter D.
%
%   Y = filter(D,X) filters input data, X, with the digital filter, D, and
%   produces output data Y. Input X can be a double or single precision
%   vector or matrix with as many columns as input channels.
%
%   fvtool(D) opens the filter visualization tool to analyze the response
%   of digital filter D.
%
%   % Example 1:
%   %   Design a highpass IIR filter with order 8, passband frequency of 
%   %   75 KHz, and a passband ripple of 0.2 dB. Sample rate is 200 KHz.
%   %   Visualize the filter response and apply the filter to a vector of
%   %   random data. 
%
%   hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%   fvtool(hpFilt) % visualize filter response
%
%   data = randn(1000,1);
%   y = filter(hpFilt,data); % apply filter to your data
%
%   % Example 2:
%   %   Design a lowpass FIR filter with normalized passband frequency of
%   %   0.25*pi rad/s, stopband frequency of 0.35*pi rad/s, passband ripple
%   %   of 0.5 dB, and stopband attenuation of 65 dB. Use a Kaiser window 
%   %   design method. Apply the filter to a vector of random data. 
%   
%   lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 0.25,...
%            'StopbandFrequency', 0.35, 'PassbandRipple', 0.5, ...
%            'StopbandAttenuation', 65, 'DesignMethod', 'kaiserwin');
%
%   data = randn(1000,1);
%   y = filter(lpFilt,data); % apply filter to your data
%
%   % Example 3:
%   %   Design a bandpass FIR filter with a lower cutoff frequency of 500
%   %   Hz and a higher cutoff frequency of 560 Hz. The sample rate is 1500
%   %   Hz. Visualize the filter response using freqz. 
%   
%   bpFilt = designfilt('bandpassfir', 'FilterOrder', 20, ...
%            'CutoffFrequency1', 500, 'CutoffFrequency2', 560,...
%            'SampleRate', 1500);
% 
%   freqz(bpFilt)
%
%   % Example 4:
%   %   Launch the Filter Design Assistant so that you can edit the filter 
%   %   designed in the example above.
%
%   designfilt(bpFilt)
%
%   See also <a href="matlab:help digitalFilter/filter">filter</a>, <a href="matlab:help digitalFilter">digitalFilter</a>.

%#ok<*NASGU>
%#ok<*AGROW>

nargoutchk(0,1)

% Detect if the function is in test mode
testModeIdx = strcmp(varargin,'TestMode');
testMode = any(testModeIdx);
varargin(strcmp(varargin,'TestMode')) = [];
doNotOfferAssistantFlag = false;
outputVarNameParseFail = false;

requestedResponse = [];
throwCommandLineError = false;

% Get caller information
fromPcodedFile = false; % code coming from pcoded file
stk =  dbstack('-completenames');
[fromCommandLine,fromAnonymousFunction] = getCallerInfo(stk,testMode);
if ~fromCommandLine
  callingFile = stk(2).file;
  callingLine = stk(2).line;  
  [~,~,callingFileExt] = fileparts(callingFile);
  fromPcodedFile = strcmp(callingFileExt,'.p');
end            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         No inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
  % No inputs, launch response selector
  
  % Verify if we have correct environment to open assistant
  checkIfAssistantSupported;
  
  varargout{1} = [];
    
  err = ''; % define empty error in case user cancels the assistant
  
  % Check if caller is a function or the command line. Sniff the output
  % variable name.
  [~,~,callingCode] = getCallerInfo(stk,testMode);
  [outputVarName,outputVarNameParseFail] = getOutputVarName(stk,callingCode,testMode);
  
  % Do not launch assistant unless we were able to parse the code  
  if outputVarNameParseFail
    error(message('signal:designfilt:CannotParseCodeToLaunchAssistant'))
  else
    % hdes is the handle to the filter desing assistant design object
    hdes = filterbuilder('FromFilterDesigner','OutputVarName',outputVarName);
    if isempty(hdes)
      % hdes will be empty if user cancels the response selector.
      return
    end  
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        digitalFilter INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif isa(varargin{1},'digitalFilter')
  
  % Verify if we have correct environment to open assistant
  checkIfAssistantSupported;
  
  if nargin > 1
    error(message('signal:designfilt:TooManyInputsWhenEditingAFilter'))
  end
  if nargout > 0
    error(message('signal:designfilt:TooManyOutputsWhenEditingAFilter'))
  end
      
  % Parse input var name using MTREE
  getParamsFromInputNameFcn = false;
  inputParamValueName = '';
  [~,~,callingCode] = getCallerInfo(stk,testMode);
  try
    [inputParamValueName,doNotOfferAssistantFlag] = getInputParamNames(stk,[],callingCode,testMode);    
    inputParamValueName = inputParamValueName{:};
  catch
    inputParamValueName = {};
    getParamsFromInputNameFcn = true;
    doNotOfferAssistantFlag = true;
  end

  if ~doNotOfferAssistantFlag && (isempty(inputParamValueName) || getParamsFromInputNameFcn)
    inputParamValueName = inputname(1);   
  end
  
  % Do not launch assistant unless we were able to parse the code
  if isempty(inputParamValueName)
    error(message('signal:designfilt:CannotParseCodeToLaunchAssistant'))
  else
    % Call filterbuilder, offer var name editor if we could not parse the
    % input variable name.
    D = varargin{1};
    hdes = filterbuilder('FromFilterDesigner', todfilt(D),...
      'Response',[D.FrequencyResponse D.ImpulseResponse],...
      'OutputVarName',inputParamValueName,'EditDigitalFilter');
    
    waitfor(hdes,'DialogClosed',true)
    
    if ~isempty(hdes.OutputCode)
      % Paste and evaluate code
      com.mathworks.mlservices.MLExecuteServices.executeCommand(hdes.OutputCode);      
    end
  end
  return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     One or more parameter inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
  
  varargout{1} = [];

  % Check valid inputs and generate fdesign code that can be used by
  % filterbuilder
  
  % Get var names from inputs
  inputParamValueNames = cell(1,numel(varargin));
  for idx = 1:numel(varargin)
    n = inputname(idx);
    inputParamValueNames{idx} = n;
  end 
  
  [err,requestedResponse,parseParams,h] = parseAndDesignFilter(inputParamValueNames, varargin{:});     
  
  % If there was an error: 
  %
  % 1) If fdesign code was generated, evaluate it to
  % make sure that we error out in cases where we have valid specification
  % sets, but with an invalid specified value (this can happen when the
  % error is due to incorrect design methods, design options, or due to
  % unknown parameters, here the generated fdesign code will contain user
  % defined parameter values). 
  %
  % 2)Parse input parameters.
  if ~isempty(err)
    if isDesktopAndJavaAndDDGOn
      
      [~,~,callingCode] = getCallerInfo(stk,testMode);
      
      % Check generated fdesign code, make sure we catch errors due to
      % invalid parameter values right here.
      if isfield(parseParams,'FdesCode') && ~isempty(parseParams.FdesCode)
        try
          fd = eval(parseParams.FdesCode);
          if isfield(parseParams,'DesignOptsStruct') && ~isempty(parseParams.DesignOptsStruct)
            tmpFilt = design(fd,parseParams.FdesMethods,parseParams.DesignOptsStruct);
          else
            tmpFilt = design(fd,parseParams.FdesMethods);
          end
        catch ME
          throwAsCaller(ME);
        end
      end
      
      % Parse input parameters
      try
        [inputParamValueNames,doNotOfferAssistantFlag] = getInputParamNames(stk,parseParams,callingCode,testMode);
      catch
        % If we cannot parse the code, then do not assume knowledge of input
        % param values. Don not offer click here links to correct code.
        inputParamValueNames = {};
        doNotOfferAssistantFlag = true;
      end
    else
      callingCode = '';
      doNotOfferAssistantFlag = true;
    end
  end
  
  if isempty(err) && ~parseParams.ResponseOnly
    % No errors and there were more inputs than just the response, so
    % return the designed filter if an output was requested.   
    varargout{1} = h;
    return;
  elseif isempty(err) && parseParams.ResponseOnly
    % No error but only the response was passed. This means we want to
    % launch the assistant without any popup dialog.
    
    % Verify if we have correct environment to open assistant
    checkIfAssistantSupported;
        
    % Check if caller is a function or the command line. Sniff the output
    % variable name.   
    [~,~,callingCode] = getCallerInfo(stk,testMode);
    [outputVarName,outputVarNameParseFail] = getOutputVarName(stk,callingCode,testMode);
    
    if doNotOfferAssistantFlag || outputVarNameParseFail
      error(message('signal:designfilt:CannotParseCodeToLaunchAssistant'))
    else
      hdes = filterbuilder('FromFilterDesigner','Response',requestedResponse,'OutputVarName',outputVarName);
    end
    
  else    
    % If there was an error then bring up the assistant and/or error out
    
    % Get output variable name from the command line or file name and line
    % if caller is a file or script (calling call is obtained above).
    [outputVarName,outputVarNameParseFail] = getOutputVarName(stk,callingCode,testMode);    
    if ~fromCommandLine
      callingFile = stk(2).file;
      callingLine = stk(2).line;
      
      [~,~,callingFileExt] = fileparts(callingFile);
      fromPcodedFile = strcmp(callingFileExt,'.p');       
    end
    
    % Check preference and if we are comming from the command line
    if ~fromCommandLine || doNotOfferAssistantFlag || ...
        outputVarNameParseFail || isdeployed
      % Never show a pop-up for the above cases Cases: From Script/Fn,
      % no-desktop, no-Java, P-coded file, compiled file, invalid code
      % parsing. Checks for these are done in getInputParamNames and
      % getOutputVarName
      launchDlg = false;
    else
      prefs = getpref('dontshowmeagain');
      if isfield(prefs,'filterDesignAssistant')
        launchDlg = ~getpref('dontshowmeagain','filterDesignAssistant');
      else
        launchDlg = true;
      end
    end
    
    % Launch popup dialog if turned-on in preferences and show assistant.
    % Hold code execution until assistant is closed. If assistant help is
    % declined, then error out on the command line.
    
    % Launch dialog if preference has not been turned off
    if launchDlg
      str1 = getString(message('signal:designfilt:FilterDesignAssistanceIsAvailable'));      
      str2 = sprintf('%s\n\n%s\n\n%s',...
        getString(message('signal:designfilt:ErrorDetectedCallingDesignFilt')),err,...
        getString(message('signal:designfilt:UseAssistantToGenerateMatlabCode')));
      
      str3 = getString(message('signal:designfilt:WouldLikeToLaunchAssistant'));
      
      choice = signal.internal.designfiltDlg({str1,str2,str3},getString(message('signal:designfilt:DoNotShowThisMessageAgain')));
      
      switch lower(choice)
        
        case 'yes'
          if isempty(requestedResponse)
            hdes = filterbuilder('FromFilterDesigner','OutputVarName',outputVarName);
            
          elseif parseParams.ResponseOnly
            
            hdes = filterbuilder('FromFilterDesigner','Response',requestedResponse,'OutputVarName',outputVarName);
            
          elseif isfield(parseParams,'FdesCode') && ~isempty(parseParams.FdesCode)
            
            codestr = sprintf('design(%s,''%s'')',parseParams.FdesCode,parseParams.FdesMethods);
            
            % Evaluate code to get a dfilt to populate the GUI with the
            % right magnitude and frequency constraint widgets.
            w = warning('off'); %#ok<WNOFF>
            tmpFilt = eval(codestr);
            warning(w)
            
            % Pass user inputs to populate GUI
            [propNames,propValues,propValueNames] = getNames(parseParams,inputParamValueNames);
            
            hdes = filterbuilder('FromFilterDesigner',tmpFilt,'Response',requestedResponse,'OutputVarName',outputVarName,...
              'PropertyNames',propNames,'PropertyValues',convertCellContentsToStrings(propValues),'PropertyValueNames',propValueNames);
            
          else
            hdes = filterbuilder('FromFilterDesigner','Response',requestedResponse,'OutputVarName',outputVarName);
          end
          
        case 'check'
          %do not show this again
          setpref('dontshowmeagain','filterDesignAssistant',true)
          throwCommandLineError = true;
          
        case {'no',[]}
          throwCommandLineError = true;
      end
    else
      throwCommandLineError = true;
    end
  end
end

% Add code to command line if needed
if ~throwCommandLineError
  
  % Pause until assistant dialog is closed. At close time the dialog will
  % set its 'DialogClosed' property to true to release the pause.
  waitfor(hdes,'DialogClosed',true)
  
  % Code will be empty if user cancels the assistant, if this happens,
  % throw command line error unless error is empty which means user just
  % passed in the response name. Otherwise, put and execute code on the
  % command line; the designed filter is held by the DfiltDesign property
  % of the designer object hdes.
  if ~isempty(hdes) && ~isempty(hdes.OutputCode)
    varargout{1} = digitalFilter(hdes.DfiltDesign);
    com.mathworks.mlservices.MLExecuteServices.executeCommand(hdes.OutputCode);
    return;
  else
    throwCommandLineError = true;
    if nargout > 0
      varargout{1} = [];
    end
  end
end

% Throw errors on the command line
if throwCommandLineError && ~isempty(err)
  if isfield(parseParams,'InputRecommendation') && ~isempty(parseParams.InputRecommendation)
    err = parseParams.InputRecommendation;
  end
  % Command line error
  if fromCommandLine
    launchStr = getString(message('signal:designfilt:AssistantToDesignFilt'));
  else
    launchStr = getString(message('signal:designfilt:AssistantThatCanCorrectYourCode'));
  end
  
  if doNotOfferAssistantFlag || outputVarNameParseFail || isdeployed 
  % Do not add click here link in the above cases, just throw the error
    error(parseParams.ErrorId,err);    
  elseif fromCommandLine
    % When user clicks on the 'click here' link and completes the design we
    % add correct line of code to the command line with arbitrary var
    % name
    outputVarName = '';
    if isempty(requestedResponse)
      error(parseParams.ErrorId,...
      	'\n<a href="matlab:filterbuilder(''FromFilterDesigner'',''AddCodeToCommandLine'',''OutputVarName'',''%s'')">Click here</a> %s.\n\n%s',...
        outputVarName,launchStr,err);
      
    elseif isfield(parseParams,'FdesCode') && ~isempty(parseParams.FdesCode)
      codestr = sprintf('design(%s,''%s'')',parseParams.FdesCode,parseParams.FdesMethods);
      
      [propNames,propValues,propValueNames] = getNames(parseParams,inputParamValueNames);
      
      error(parseParams.ErrorId,...
        '\n<a href="matlab:filterbuilder(''FromFilterDesigner'',''AddCodeToCommandLine'',''Response'',''%s'',''OutputVarName'',''%s'',''PropertyNames'',{%s},''PropertyValues'',{%s},''PropertyValueNames'',{%s},%s)">Click here</a> %s.\n\n%s',...
        requestedResponse,outputVarName,cell2str(propNames),cell2str(propValues),cell2str(propValueNames),codestr,launchStr,err);
    else
      error(parseParams.ErrorId,...
        '\n<a href="matlab:filterbuilder(''FromFilterDesigner'',''AddCodeToCommandLine'',''Response'',''%s'',''OutputVarName'',''%s'')">Click here</a> %s.\n\n%s',...
        requestedResponse,outputVarName,launchStr,err);
    end
    
  else
    % From File/Script
    % When user clicks on the 'click here' link and completes the design we
    % open the file and correct the erroneous line of code
    
    if fromAnonymousFunction
    % Do not offer code correction if anonymous as we do not know what line 
    % it is coming from  
       error(parseParams.ErrorId,err);      
    end
    
    if isempty(requestedResponse)
       error(parseParams.ErrorId,...
        '\n<a href="matlab:filterbuilder(''FromFilterDesigner'',''CorrectCode'',''FileName'',''%s'',''FileLine'',%s)">Click here</a> %s\n\n.%s',...
        callingFile,mat2str(callingLine),launchStr,err);

    elseif isfield(parseParams,'FdesCode') && ~isempty(parseParams.FdesCode)
      codestr = sprintf('design(%s,''%s'')',parseParams.FdesCode,parseParams.FdesMethods);
        
      [propNames,propValues,propValueNames] = getNames(parseParams,inputParamValueNames);
      
      error(parseParams.ErrorId,...
        '\n<a href="matlab:filterbuilder(''FromFilterDesigner'',''Response'',''%s'',''CorrectCode'',''FileName'',''%s'',''FileLine'',%s,''PropertyNames'',{%s},''PropertyValues'',{%s},''PropertyValueNames'',{%s},%s)">Click here</a> %s.\n\n%s',...
        requestedResponse,callingFile,mat2str(callingLine),cell2str(propNames),cell2str(propValues),cell2str(propValueNames),codestr,launchStr,err);
        
    else
      error(parseParams.ErrorId,...
        '\n<a href="matlab:filterbuilder(''FromFilterDesigner'',''Response'',''%s'',''CorrectCode'',''FileName'',''%s'',''FileLine'',%s)">Click here</a> %s.\n\n%s',...
        requestedResponse,callingFile,mat2str(callingLine),launchStr,err);
    end
  end
end

%--------------------------------------------------------------------------
function [err,requestedResponse,parseParams,h] = parseAndDesignFilter(inputValueNames, varargin)
%checkfilter Parse inputs and design filter if inputs are correct

% inputValueNames contains names used by user instead of numeric values
% (for example Fs to denote the value of sample rate)

if numel(inputValueNames) ~= numel(varargin)
  error(message('signal:designfilt:LengthOfInputMustEqualLengthOfVarargin'));
end

% Default values for input args:
propNames = {};
propValues = {};
if isempty(varargin)
  filterType = '';
else
  filterType = varargin{1};      
  if nargin > 1
    if mod(length(varargin(2:end)),2)~= 0
      error(message('signal:designfilt:InputsMustBePVP'));
    end    
    propNames = varargin(2:2:end);
    
    checkForRepeatedProperties(propNames);
    
    if ~iscellstr(propNames)
      error(message('signal:designfilt:PropertyMustBeStrings'));
    end    
    propValues = varargin(3:2:end);
    inputValueNames = inputValueNames(3:2:end);
  end
end

% Check that we have correct set of specifications
parserObj = signal.internal.DesignfiltProcessCheck;
[err,requestedResponse,parseParams] = checkConstraints(parserObj,filterType,propNames,propValues,inputValueNames);

if ~isempty(err) || parseParams.ResponseOnly 
  h = [];
  return;
end
  
% Design filter
h = parseParams.DesignedFilter;

%--------------------------------------------------------------------------
function [propNames,propValues,propValueNames] = getNames(parseParams,inputParamValueNames)

propNames = {};
propValues = {};
propValueNames = {};

propNames = [propNames parseParams.PropertyNames];
propNames = [propNames parseParams.SampleRateProps];
propNames = [propNames parseParams.DesignMethodProps];
propNames = [propNames parseParams.DesignOptsProps];

propValues = [propValues parseParams.PropertyValues];
propValues = [propValues parseParams.SampleRateVals];
propValues = [propValues parseParams.DesignMethodVals];
propValues = [propValues parseParams.DesignOptsVals];

propValueNames = [propValueNames parseParams.PropertyValNames];
propValueNames = [propValueNames parseParams.SampleRateValNames];
propValueNames = [propValueNames parseParams.DesignMethodValNames];
propValueNames = [propValueNames parseParams.DesignOptsValNames];

% Replace propValueNames with ones we parsed from the code. 
if ~isempty(inputParamValueNames)
  for k = 1:numel(propNames)
    idx = find(strcmp(inputParamValueNames,propNames{k}));
    if ~isempty(idx)
      propValueNames{k} = inputParamValueNames{idx+1};
    end
  end
end

%--------------------------------------------------------------------------
function [s,N] = cell2str(c,d)
%cell2str Return delimiter-separated string.

if nargin<2
  d = ',';
end
N = numel(c);
s = '';
for idx=1:N
  if ischar(c{idx})
    if strncmp(c{idx},'{',1) && isempty(strfind(c{idx},'@'))
      s = [s '''' strrep(c{idx},'''','''''') ''''];
    else
      s = [s '''' c{idx} '''']; 
    end
  elseif islogical(c{idx})
    s = [s mat2str(c{idx})];  
  elseif isnumeric(c{idx})
    s = [s '''' mat2str(c{idx}) '''']; 
  elseif isa(c{idx},'function_handle')
    s = [s '@' func2str(c{idx})];
  elseif iscell(c{idx})
    % Must be window cell
    c{idx} = convertCellToString(c{idx});
    s = [s c{idx}];
  end
  if idx<N
    s = [s d]; 
  end
end
%--------------------------------------------------------------------------
function c = convertCellContentsToStrings(c)
%convertCellContentsToStrings Convert cell contents to strings

for idx = 1:numel(c)
  if isnumeric(c{idx})
    c{idx} = mat2str(c{idx});
  elseif islogical(c{idx})
    c{idx} = mat2str(c{idx});
  elseif isa(c{idx},'function_handle')
    c{idx} = ['@' func2str(c{idx})];
  elseif ischar(c{idx})
    c{idx} = ['''' c{idx} ''''];     
  elseif iscell(c{idx})
    % Must be window cell
    c{idx} = convertCellToString(c{idx});   
  end
end
%--------------------------------------------------------------------------
function  c = convertCellToString(c)
%convertCellToString Convert a cell to a string '{'a',b}'
newC = c;
k = cellfun('isclass',c,'function_handle');
if any(k)
  newC{k} = ['@' func2str(c{k})];
end
k = cellfun('isclass',c,'char');
if any(k)
  newC{k} = ['''' c{k} ''''];
end
str = '{';
for k = 1:numel(newC)
  if isnumeric(newC{k})
    newC{k} = mat2str(newC{k});
  end
  str = [str newC{k}];
  if k < numel(newC)
    str = [str ','];
  end
end
str = [str '}'];
c = str;
%--------------------------------------------------------------------------
function [outputVarName,errFlag] = getOutputVarName(stk,callingCode,testMode)
%getOutputVarName This function returns the output var name if caller comes
%from command line. Parsing of code from a file is done later in
%filterbuilder. 

outputVarName = '';
errFlag = false;

[fromCommandLine,fromAnonymousFunction] = getCallerInfo(stk,testMode);

if ~isDesktopAndJavaAndDDGOn() || fromAnonymousFunction
  % Do not parse code if we have no desktop or if we are called from an
  % anonymous function.
  errFlag = true;
  return;
  
elseif fromCommandLine
  % If call comes from command line and we have a desktop, then explore the
  % command line text, or the command line history and extract the output
  % variable name. If call comes from a function/script, the code will be
  % analyzed in filterbuilder.m
  
  designfiltCallProps = signal.internal.parseDesignfiltCallCodeFromCmdLine(callingCode);
  
  if designfiltCallProps.ErrorFlag
    errFlag = true;
    return;
  end
  
  outputVarName = designfiltCallProps.OutputVarName;
end
%--------------------------------------------------------------------------
function [inputParamValueNames,doNotOfferAssistantFlag] = getInputParamNames(stk,parseParams,callingCode,testMode)
% Input stk is the dbstack

doNotOfferAssistantFlag = false;
inputParamValueNames = {};

fromAnonymousFunction = false;
fromCommandLine = (numel(stk)==1) || (numel(stk)==2 && isempty(stk(2).file));

if numel(stk)>=2 && isempty(stk(2).file)
  fromAnonymousFunction = true;
end

if testMode
  fromCommandLine = true;
end

if ~isDesktopAndJavaAndDDGOn() || fromAnonymousFunction
  % Do not parse code if we have no desktop or if we are called from an
  % anonymous function.
  doNotOfferAssistantFlag = true;
  return;
elseif fromCommandLine    
  % If call comes from command line then explore the desktop text. If we
  % cannot get the desktop text, then explore the command line history.  
      
  designfiltCallProps = ...
    signal.internal.parseDesignfiltCallCodeFromCmdLine(callingCode,'all'); 
else
  % Call comes from a script/function
  callingFile = stk(2).file;
  callingLine = stk(2).line;
  
  [~,~,callingFileExt] = fileparts(callingFile);
  fromPcodedFile = strcmp(callingFileExt,'.p');
  
  if fromPcodedFile
    doNotOfferAssistantFlag = true;
    return;
  end
  
  designfiltCallProps  = signal.internal.parseDesignfiltCallCodeFromFile(...
    callingFile,callingLine);  
end

if designfiltCallProps.ErrorFlag 
  % Could not parse code, so return with empty set of param value names
  if ~fromCommandLine
    doNotOfferAssistantFlag = true;    
  end
  return;
end

% Get input param names and remove ambiguities
inputParamNames = designfiltCallProps.InputParamNames;
inputParamNames(strcmp(inputParamNames,'''TestMode''')) = [];
inputParamValues = designfiltCallProps.InputParamValues;
inputParamValuesAreStrs = designfiltCallProps.InputParamValuesAreStrings;

% Concatenate values. Remove response. Remove double quotes from params
% that are truly input strings. 
inputParamValueNames = strrep(inputParamValues(1),'''',''); % response
for idx = 2:numel(inputParamNames)
  currentParamValue = inputParamValues{idx};
  if inputParamValuesAreStrs{idx}
    currentParamValue = strrep(currentParamValue,'''','');
  end
  
  name = strrep(inputParamNames{idx},'''','');
  if ~isempty(parseParams)
    % Need to call this with more than 3 outputs, otherwise, function
    % errors out internally if it finds an abiguity or an invalid
    % parameter. This function checks for matches with a list, and if a
    % match is found without ambiguity it sets the output to the list
    % value. 
    [name,~,~,~,~] = signal.internal.DesignfiltProcessCheck.validateStringList(...
      parseParams.FullInputNameList,name);
  end  
  % Correct strings for windows, design methods, and match exactly
  if strcmpi(name,'Window') 
    list = signal.internal.DesignfiltProcessCheck.getWindowsList;
    if inputParamValuesAreStrs{idx}            
      [currentParamValue,~,~,~,~] = signal.internal.DesignfiltProcessCheck.validateStringList(...
        list,currentParamValue);
      currentParamValue = currentParamValue{:};
    elseif strncmp(currentParamValue,'{',1)
      kidx = strfind(currentParamValue,'''');
      if numel(kidx) == 2 && isempty(strfind(currentParamValue,'@'))
        [winName,~,~,~,~] = signal.internal.DesignfiltProcessCheck.validateStringList(...
        list,currentParamValue(kidx(1)+1:kidx(2)-1));
        currentParamValue = [currentParamValue(1:kidx(1)) winName{:} currentParamValue(kidx(2):end)];      
      end      
    end                 
  elseif strcmpi(name,'DesignMethod')
    list = signal.internal.DesignfiltProcessCheck.getDesignMethodsList;
    [currentParamValue,~,~,~,~] = signal.internal.DesignfiltProcessCheck.validateStringList(...
      list,currentParamValue);
    currentParamValue = currentParamValue{:};
  elseif strcmpi(name,'MatchExactly')
    [currentParamValue,~,~,~,~] = signal.internal.DesignfiltProcessCheck.validateStringList(...
      {'passband','stopband','both'},currentParamValue);    
    currentParamValue = currentParamValue{:};
  end
  
  inputParamValueNames = [inputParamValueNames,name{:}, currentParamValue];
end
%--------------------------------------------------------------------------
function [fromCommandLine,fromAnonymousFunction, callingCode] = getCallerInfo(stk,testMode)
%getCallerInfo 
callingCode = '';
fromAnonymousFunction = false;
fromCommandLine = (numel(stk)==1) || (numel(stk)==2 && isempty(stk(2).file));

if testMode
  fromCommandLine = true;
end

if numel(stk)>=2 && isempty(stk(2).file)
  fromAnonymousFunction = true;
end

if fromCommandLine && nargout > 2 && isDesktopAndJavaAndDDGOn()
  callingCode = getCallHistory(); 
end

%--------------------------------------------------------------------------
function checkForRepeatedProperties(propNames)
%checkForRepeatedProperties Error out if properties have been specified
%more than once

repeatedPropList = {};
propNamesUnique = unique(propNames,'stable');
for idx = 1:numel(propNamesUnique)
  if sum(ismember(propNames,propNamesUnique{idx})) > 1
    repeatedPropList = [repeatedPropList propNamesUnique{idx}];
  end
end

if ~isempty(repeatedPropList)  
  N = numel(repeatedPropList);
  s = '';
   d = ', ';
  for i=1:N
    s = [s repeatedPropList{i}];
    if i<N
      if i == (N-1) && N>1
        s = [s d 'and '];
      else
        s = [s d];
      end
    end
  end
  
  if N > 1
    error(message('signal:designfilt:RepeatedParams',s));
  else
    error(message('signal:designfilt:RepeatedParam',s));
  end  
end
%--------------------------------------------------------------------------
function callHist = getCallHistory()
% Get call from command line or command line history
try
  jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
  cmdWin = jDesktop.getClient('Command Window');
  jTextArea = cmdWin.getComponent(0).getViewport.getComponent(0);
  callHist = char(get(jTextArea,'Text'));  
catch
  callHist = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
end

%--------------------------------------------------------------------------
function checkIfAssistantSupported()
% checkIfAssistantSupported Verify if we can offer the assistant
% User must have desktop, Java, and not have called designfilt from a
% compiled file. 

if ~isDesktopAndJavaOn
  % Assistant needs desktop and Java.
  error(message('signal:designfilt:AssistantRequiresDesktopAndJava'))
elseif isdeployed
  % Assistant not supported if script is compiled
  error(message('signal:designfilt:AssistantNotSupportedInCompiled'));
elseif ~isDesktopAndJavaAndDDGOn
  % Assistant needs DDG support
  error(message('signal:designfilt:AssistantNotSupportedInMOTW'));
end
%--------------------------------------------------------------------------
function flag = isDesktopAndJavaOn()
%isDesktopAndJavaOn Check if desktop and java are available

isDesktopOn = usejava('desktop');
isJavaOn = all([usejava('jvm') usejava('swing') usejava('awt')]);
flag = isDesktopOn && isJavaOn;
 
%--------------------------------------------------------------------------
function flag = isDesktopAndJavaAndDDGOn()
%isDesktopAndJavaAndDDGOn Check if desktop and java are available and if
%DDG is supported

isDDGOn = signal.internal.SPTCustomSettings.isDDGSupported;    
flag = isDesktopAndJavaOn && isDDGOn;


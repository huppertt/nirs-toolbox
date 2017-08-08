function [outputVarName, errFlag, codeCell, fcallLineIdx, indentLength,inputParamInfo] = parseDesignfiltCallCode(cmdSet)
%parseDesignfiltCallCode Find output variable name and parse code
%involved in designfilt call

% cmdSet can be a cell array with code lines, a java.lang.String[]
% object, or a string.

%   Copyright 2013 The MathWorks, Inc.

% Get sting of code involving the designfilt call
if ischar(cmdSet)
  str = cmdSet;
else
  str = getString(cmdSet);
end

% Parse the code
if nargout > 2
  [outputVarName, errFlag, codeCell, fcallLineIdx, indentLength, inputParamNames, inputParamValues, inputParamValuesAreStrs] = parseCode(str);
  inputParamInfo.Names = inputParamNames;
  inputParamInfo.Values = inputParamValues;
  inputParamInfo.ValuesAreStrs = inputParamValuesAreStrs;
else
  [outputVarName, errFlag] = parseCode(str);
end

%--------------------------------------------------------------------------
function [str,strCell] = getString(cmdSet)

for idx = 1:numel(cmdSet)
  
  if isa(cmdSet,'java.lang.String[]') || isa(cmdSet,'java.lang.String')
    strCell{idx} = toCharArray(cmdSet(idx)).'; %#ok<*AGROW>
  else
    strCell{idx} = cmdSet{idx};
  end
  
  % Concatenate a string, do not use char since we need a one dimensional
  % string to pass to mtree
  if idx == 1
    str = strCell{idx};
  else
    str = sprintf('%s\n%s',str,strCell{idx});
  end
end

%--------------------------------------------------------------------------
function [outputVariableName, errFlag, codeCell, fcallLineIdx, indentLength, inputParamNames, inputParamValues, inputParamValuesAreStrs] = parseCode(str)
%parseCode Parse output variable name and code involving the designfilt
%call using mtree

outputVariableName = '';
codeCell = {};
fcallLineIdx = [];
indentLength = [];
errFlag = false;
inputParamNames = {};
inputParamValues = {};
inputParamValuesAreStrs = {};

T = mtree(str);

if isempty(T)
  errFlag = true;
  return;
end

if ~isempty(mtfind(T, 'Kind', 'ERR'))
  errFlag = 'incompleteCode';
  return;
end

P = mtfind(T, 'Fun','designfilt');

if isempty(P) || count(P) > 1 || ~isempty(mtfind(full(P),'Kind','ANON'))
  % If we did not find a designfilt call, or we find more than one call, or
  % if call comes from anonymous function, then return with error flag set
  % to true.
  errFlag = true;
  return;
end

outputVariableName = '';
done = isempty(P);
while ~done
  P = P.trueparent;
  
  if ~isempty(P)
    E = mtfind( P, 'Kind', 'EQUALS');
    C = mtfind( P, 'Kind', 'CALL');
    if ~isempty(C)
      designfiltTree = Full(C);
    end
    if ~isempty(E)
      outputVariableName = convertTree2str(E.Left);
      designfiltTree = Full(E);
      done = true;
    end
    
  else
    done = true;
  end
end

if nargout > 2
  fcall = tree2str(designfiltTree); % parse designfilt call line
  newStr = tree2str(T); % parsed code from the input string
  codeCell = regexp(newStr, '[\f\n\r]', 'split');
  
  for idx = 1:numel(codeCell)
    flag = strfind(codeCell{idx},fcall);
    if ~isempty(flag)
      fcallLineIdx = idx;
      indentLength = flag-1;
    end
  end
  
  % Output vat name. Go through all the parameters to the right of the
  % designfilt call.
  P = full(mtfind(T, 'Fun','designfilt'));
  
  X = P.Arg.Right.Right; % first input param (response)
  if isempty(X) || ~strcmpi(kinds(trueparent(X)),'call')
    % No equal sign
    X = P.Arg.Right;
  end
  if isempty(X.Next) && ~strcmpi(kinds(X),'string')
    % Only ine input, we are trying to edit a digital filter    
    inputParamNames = {''};
    inputParamValues = {convertTree2str(X)};
    inputParamValuesAreStrs = {false};
    return;
  else
    % We are trying to edit a digital filter
    inputParamNames = {'FilterResponse'};
    inputParamValues = {convertTree2str(X)};
    inputParamValuesAreStrs = {true};   
  end
  
  isDone = false;
  cnt = 1;
  while ~isDone
    X = X.Next;
    if isempty(X)
      isDone = true;
    else
      if mod(cnt,2) == 1
        if strcmpi(kinds(X),'string')
          inputParamNames = [inputParamNames; {convertTree2str(X)}];
        else
          inputParamNames = [inputParamNames; {' '}];
        end
      else
        inputParamValues = [inputParamValues; {convertTree2str(X)}];
        inputParamValuesAreStrs = [inputParamValuesAreStrs {strcmpi(kinds(X),'string')}];
      end      
      cnt = cnt + 1;      
    end
  end  
end

%--------------------------------------------------------------------------
function str = convertTree2str(T)
%convertTree2str Get string from tree T and remove blank spaces
str = tree2str(T);
str = strrep(str,' ','');


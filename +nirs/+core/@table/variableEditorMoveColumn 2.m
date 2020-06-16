function [moveCode,msg] = variableEditorMoveColumn(this,varName,startCol,endCol)
%   This function is undocumented and will change in a future release


%   Copyright 2011-2012 The MathWorks, Inc.

msg = '';
[~,varIndices] = variableEditorColumnNames(this);
startIndex = find(varIndices<=startCol,1,'last');
endIndex = find(varIndices<=endCol,1,'last');

if startIndex<endIndex
    if startIndex>=3 % x = [x(:,[1:startIndex-1
        moveCode = [varName ' = ' varName '(:,[1:' num2str(startIndex-1) ' '];
    elseif startIndex==2 % x = [x(:,[1 
        moveCode = [varName ' = ' varName '(:,[1 '];
    else % x = x(:,[ 
        moveCode = [varName ' = ' varName '(:,['];
    end
    if endIndex-1>startIndex+1 % ... startIndex+1:endIndex-1) ...
        moveCode = [moveCode ...
                    num2str(startIndex+1) ':' num2str(endIndex-1) ' '];
    elseif endIndex-1==startIndex+1 % ... startIndex+1 ...
        moveCode = [moveCode ... 
                    num2str(startIndex+1) ' '];
    end
    if endIndex==size(this,2) % ... startIndex end])
        moveCode = [moveCode ...
                    num2str(startIndex) ' end]);']; 
    elseif endIndex==size(this,2)+1 % ... startIndex])
        moveCode = [moveCode ...
                    num2str(startIndex) ']);']; 
    else % ... startIndex endIndex:end])
        moveCode = [moveCode ...
                    num2str(startIndex) ' '  num2str(endIndex) ...
                   ':end]);']; 
    end
elseif endIndex<startIndex
    if endIndex>=3 % x = x(:,1:endIndex-1 
        moveCode = [varName ' = ' varName '(:,[1:' num2str(endIndex-1) ' '];
    elseif endIndex==2 % x = [x(:,1)  
        moveCode = [varName ' = ' varName '(:,[1 '];
    else % x = x(:,[ 
        moveCode = [varName ' = ' varName '(:,['];
    end
    
    % ... startIndex ...
    moveCode = [moveCode ... 
                num2str(startIndex) ' ']; 
                
    if startIndex-1>endIndex  % ... endIndex:startIndex-1 ...          
        moveCode = [moveCode ...
                    num2str(endIndex) ':' num2str(startIndex-1)];
    else % ... endIndex ... 
         moveCode = [moveCode ...
                     num2str(endIndex)];
    end
    
    %
    if startIndex+1<size(this,2)
        moveCode = [moveCode ...
                    ' ' num2str(startIndex+1) ':end]);'];
    elseif startIndex+1==size(this,2)            
        moveCode = [moveCode ...
                    ' end]);'];  
    else
        moveCode = [moveCode ...
                    ']);'];
    end
end
    

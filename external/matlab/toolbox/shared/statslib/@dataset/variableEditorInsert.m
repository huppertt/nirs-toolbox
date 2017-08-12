function out = variableEditorInsert(this,orientation,row,col,data)
%   This function is undocumented and will change in a future release

% Performs a insert operation on data from the clipboard.

%   Copyright 2011-2012 The MathWorks, Inc.

numVars = size(this,2);
numRows = size(this,1);

% Get the inserted data as a dataset
if isa(data,'dataset')
    varData = data;
elseif iscell(data)        
    % Try to create a dataset to numeric data first
    try
        varData = dataset(cell2mat(data));
    catch me %#ok<NASGU>
        varData = dataset(data);
    end
else
    varData = dataset(data);
end

[~,varIndices] = variableEditorColumnNames(this);
insertIndex = find(varIndices<=col,1,'last');

% Find the insertion index
if strcmp('columns',orientation)
        
    % Add the pasted dataset after the last column
    %this(row:size(data,1)+row-1,end+1:end+size(data,2)) = data;
    newVarIndices = size(this,2)+1:size(this,2)+size(varData,2);
    this = subsasgn(this,struct('type','()','subs',{{row:size(varData,1)+row-1,...
        newVarIndices}}),varData); 
    for k=1:length(newVarIndices)
        % Try to assign the variable name of each inserted dataset

        % Append digits to guarantee a unique variable name
        newVarName = varData.Properties.VarNames{k};
        ind = 1;
        while any(strcmp(newVarName,this.Properties.VarNames))
            newVarName = [varData.Properties.VarNames{k},num2str(ind)];
            ind=ind+1;
        end

        %this.Properties.VarNames{newVarIndices(k)} = newVarName;
        varNameAssignment = struct('type',{'.','.','{}'},'subs',{'Properties','VarNames',{newVarIndices(k)}});
        this = subsasgn(this,varNameAssignment,newVarName);        
    end

    % Move the appended dataset to the inserted index
    %this = this(:,[1:index numVars+1:end index+1:numVars]);
    this = subsref(this,struct('type','()','subs',{{1:size(this,1),...
        [1:insertIndex-1 numVars+1:size(this,2) insertIndex:numVars]}}));

else
    if row>numRows+1
       % Add the pasted dataset after the last row
       % this(row:row+size(varData,1)-1:insertIndex+insertIndex+size(data,2)-1)
       % = varData
       out = subsasgn(this,struct('type','()','subs',...
           {{row:row+size(varData,1)-1,insertIndex:insertIndex+size(data,2)-1}}),...
           varData);  
       return
    end

   % Add the pasted dataset after the last row
   % this(numRows+1:numRows+size(varData,1):insertIndex+insertIndex+size(data,2)-1)
   % = varData
   this = subsasgn(this,struct('type','()','subs',...
       {{numRows+1:numRows+size(varData,1),insertIndex:insertIndex+size(data,2)-1}}),...
       varData);
   
   % Try to preserve the obsNames so that rows can be cut and inserted
   % with their ObsNames intact.
   if ~isempty(varData.Properties.ObsNames)
       obsNames = this.Properties.ObsNames;
       obsNames(numRows+1:numRows+size(varData,1)) = varData.Properties.ObsNames;
       obsNames = matlab.lang.makeUniqueStrings(obsNames,numRows:length(obsNames),namelengthmax);    
       this.obsnames = obsNames;
   end
   
   % Move the appended dataset to the inserted row
   % this = this([1:row numRows+1:end row+1:numRows]);
   if insertIndex<=numRows
       this = subsref(this,struct('type','()','subs',{{
          [1:row-1 numRows+1:numRows+size(varData,1) row:numRows],1:numVars}}));
   end


end


out = this;

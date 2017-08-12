function out = variableEditorPaste(this,rows,columns,data)
%   This function is undocumented and will change in a future release

% Performs a paste operation on data from the clipboard which was not
% obtained from another dataset array.

%   Copyright 2011-2013 The MathWorks, Inc.

if isa(data,'dataset')
    [nrows,ncols] = variableEditorGridSize(data);
else
    ncols = size(data,2);
    nrows = size(data,1);
end
[varNames,varIndices,classes] = variableEditorColumnNames(this);
isCatClass = cellfun(@(x)isa(x,'categorical'),this.data);
datasetSize = size(this);

% Find the dataset indices which correspond to the columnIntervals. If the
% number of pasted columns does not match the number of selected columns,
% just paste columns starting at the left-most column 
if length(columns)~=ncols
    I = unique(find(varIndices(1:end-1)>=columns(1) & varIndices(1:end-1)<=columns(1)+ncols-1));
    columns = columns(1):columns(1)+ncols-1;
else
    % Find the indices of each column
    colIndices=zeros(varIndices(end),1);
    colIndices(varIndices)=1;
    colIndices = cumsum(colIndices);
    colIndices = colIndices(1:end-1);
    
    % Index into it to convert (ignoring columns beyond the last dataset
    % column)
    I = unique(colIndices(columns(columns<=varIndices(end-1))));
end

% If the number of pasted rows does not match the number of selected rows,
% just paste rows starting at the top-most row 
if length(rows)~=nrows
    rows = rows(1):rows(1)+nrows-1;
end
 

% Paste data one variable at a time, converting from a cell array for
% non-cell dataset variables

% Paste data onto existing dataset variables
col = 1;
if ~isempty(I)
    % Loop through variable indices from min(I) to max(I)
    for k=1:length(I)
        if isa(data,'dataset')
            s = struct('type',{'()'},'subs',{{':',k}});      
            variableData = subsref(data,s);%data(:,k);
            % this(rows,I(k)) = variableData
            s = struct('type',{'()'},'subs',{{rows,I(k)}});
            this = subsasgn(this,s,variableData);
        else
            variableData = data(:,col:min(size(data,2),col+varIndices(I(k)+1)-varIndices(I(k))-1));
            col = col+varIndices(I(k)+1)-varIndices(I(k));

            % this.VarName(rows,1:size(variableData,2)) = dataset(variableData)
            s = struct('type',{'.','()'},'subs',{varNames{I(k)},{rows,...
                1:size(variableData,2)}});       
            if strcmp(classes{I(k)},'cell')
                this = subsasgn(this,s,variableData);
            elseif isCatClass(I(k)) % categorical or its subclasses
                % This supports assignment from a categorical or a cellstr.
                % Assignment of a categorical will fail if categories don't
                % match and the target is ordinal or protected.
                this = subsasgn(this,s,variableData);
            else
                try
                    this = subsasgn(this,s,cell2mat(variableData));
                catch %#ok<CTCH>
                    this = subsasgn(this,s,variableData);
                end
            end
        end
    end  
end

if isa(data,'dataset') && ~isempty(data.obsnames)
   obsNames = this.obsnames;
   extendRows = (rows > datasetSize(1));
   if isempty(obsNames)
       % Create default names for the destination and copy over all the source names.
       obsNames = dfltobsnames(1:max(max(rows),datasetSize(1)));
       obsNames(rows) = data.obsnames;
       obsNames = matlab.lang.makeUniqueStrings(obsNames,min(rows):length(obsNames),namelengthmax);
   elseif any(extendRows)
       % Preserve the destination's names, but paste names from the
       % source for obs that extend the destination
       obsNames(rows(extendRows)) = data.obsnames(extendRows);
       obsNames = matlab.lang.makeUniqueStrings(obsNames,datasetSize(1):length(obsNames),namelengthmax);
   end
   this.obsnames = obsNames;
end
   
% Any remaining pasted data should be added in the columns to the right of
% the dataset. Discontiguous column selection in this region must be
% ignored since the dataset does not allow empty variables.
appendedColumnCount = sum(columns>=varIndices(end));
if appendedColumnCount>0    
    %appendedCols = startColumn+ncols-varIndices(end);
    
    % Sub-reference the pasted data to get the part to be pasted after the 
    % dataset.
    if isa(data,'dataset')
        [~,srcVarIndices] = variableEditorColumnNames(data);
        firstAppendedColumnPosition = ncols-appendedColumnCount+1;
        srcStartIndex = find(srcVarIndices(1:end-1)<=firstAppendedColumnPosition,1,'last');
        appendedIndices = size(data,2)-srcStartIndex+1;
        s = struct('type',{'()'},'subs',{{':',size(data,2)-appendedIndices+1:size(data,2)}});
        data = subsref(data,s);
    else
        appendedIndices = appendedColumnCount;
        data = data(:,end-appendedColumnCount+1:end);
    end
    
    lastVariableIndex = length(varIndices)-1;
    for index=size(data,2)-appendedIndices+1:size(data,2) 
        % variableData = data(:,index);
        s = struct('type',{'()'},'subs',{{':',index}});
        variableData = subsref(data,s);
        
        % this(rows,lastVariableIndex+index) = variableData
        s = struct('type','()','subs',{{rows,lastVariableIndex+index}});
        
        if isa(data,'dataset')
            this = subsasgn(this,s,variableData);
               
            %this.varnames{end} = newVarName;
            newVarName = localGetUniqueVarName(data.varnames{index},this.varnames);
            varNameAssignment = struct('type',{'.','.','{}'},'subs',{'Properties','VarNames',{size(this,2)}});
            this = subsasgn(this,varNameAssignment,newVarName);
            
        else   
            try
                % Prevent char arrays
                if all(cellfun(@(x) isnumeric(x),variableData))
                    this = subsasgn(this,s,dataset(cell2mat(variableData)));
                else
                    this = subsasgn(this,s,dataset(variableData));
                end
            catch %#ok<CTCH>
                this = subsasgn(this,s,dataset(variableData));
            end
        end
    end
end
out = this;


function varName = localGetUniqueVarName(newVarName,varNames)

% Append digits to guarantee a unique variable name
ind = 1;
varName = newVarName;
while any(strcmp(varName,varNames))
    varName = [newVarName,num2str(ind)];
    ind=ind+1;
end

function [out,warnmsg] = variableEditorRowDeleteCode(this,varName,rowIntervals)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to delete rows in positions defined
% by the 2 column rowIntervals matrix. It is assumed that row intervals
% are disjoint, in monotonic order, and bounded by the number of rows 
% in the table.

%   Copyright 2011-2012 The MathWorks, Inc.

warnmsg = '';
if size(rowIntervals,1)==1
    rowString = localGenerateIntervalString(rowIntervals,size(this,1));
else
    rowString = '[';
    for k=1:size(rowIntervals,1)
        rowString = [rowString localGenerateIntervalString(rowIntervals(k,:),size(this,1))]; %#ok<AGROW>
        if k<size(rowIntervals,1)
            rowString = [rowString ',']; %#ok<AGROW>
        end
    end
    rowString = [rowString ']'];
end
out = [varName '(' rowString ',:) = [];'];

function intervalString = localGenerateIntervalString(intervalArray,len)

% Generate a string representation of a single interval where the number of
% rows or columns is len.
if intervalArray(1)==intervalArray(2)
    if intervalArray(2)==len
        intervalString = 'end';
    else
        intervalString = num2str(intervalArray(1));
    end
elseif intervalArray(2)==len
    intervalString = [num2str(intervalArray(1)) ':end'];
else
    intervalString = [num2str(intervalArray(1)) ':' num2str(intervalArray(2))];
end

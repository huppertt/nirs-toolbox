function varargout = variableEditorGridSize(a)
%   This function is undocumented and will change in a future release

% Undocumented method used by the Variable Editor to determine the number 
% of rows and columns needed to display the table data. The table
% is assumed to have 2 dimensions.

%   Copyright 2011-2014 The MathWorks, Inc.
   
varWidths = cellfun(@(x) size(x,2)*ismatrix(x)*~ischar(x)*~isa(x,'dataset')*~isa(x,'table')...
    +ischar(x)+isa(x,'dataset')+isa(x,'table'),a.data);
if (isempty(a))  
    gridSize = [0 0];
    if size(a,2)  > 0  || ~all(varWidths > 0)   
       gridSize = [-1 -1];         
    end
        
elseif all(varWidths>0)
    gridSize = [size(a,1) sum(varWidths)];
else % The Variable Editor should not use a 2d grid to display ND arrays
    gridSize = [size(a,1) 0];
end
if nargout==2
    varargout{1} = gridSize(1);
    varargout{2} = gridSize(2);
elseif nargout==1
    varargout{1} = gridSize;
end
   

   


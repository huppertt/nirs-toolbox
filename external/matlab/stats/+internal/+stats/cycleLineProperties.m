function varargout = cycleLineProperties(nG,varargin)
%Cycle through the line properties by expansion/shrinkage
%   varargout = CYCLELINEPROPERTIES(nG,varargin) expands/shrinks the values
%   in each line property given in varargin to the number of groups nG.
%   Used by plotting functions with a grouping variable such as
%   scatterhist, gscatter. Can be used on the following line properties:
%      Property           value
%      - 'LineWidth'      scalar or a numerical vector
%      - 'LineSytle'      string or cell string
%      - 'MarkerSize'     scalar or a numerical vector
%      - 'MarkerSymbol'   string or cell string
%      - 'Color'          string , cell string or a RGB color matrix. 
%                        
%   Example:
%   >> clr = 'rgb';
%   >> ls = {'-','-.',':'};
%   >> lw  = 2;
%   >> [clr, ls, lw] = internal.stats.cycleLineProperties(5,clr,ls,lw);
%   The results are:   
%   clr = 'rgbrg'
%   ls =  {'-'    '-.'    ':'    '-'    '-.'}
%   lw =  [2     2     2     2     2]

%   Copyright 2012 The MathWorks, Inc.

narginchk(2,inf);
numProps = nargin-1;
nargoutchk(0,numProps);
varargout = varargin;

for i = 1:numProps 
    if isrow(varargout{i})
        varargout{i} = expandCol(varargout{i},nG);
    elseif iscolumn(varargout{i})
        varargout{i} = expandCol(transpose(varargout{i}),nG);
        varargout{i} = transpose(varargout{i});
    elseif ismatrix(varargout{i}) % Expand rows of a color matrix
        clrMat = varargout{i};
        clrMat = expandCol(clrMat',nG);
        varargout{i} = clrMat';
    end   
end

end

function B = expandCol(A,total)
%EXPANDMAT Expand the columns of an array by replicating

colA = size(A,2);
m = ceil(total/colA);
rowB = m*colA;

B = repmat(A,1,m);

if total <rowB
    B(:,total+1:rowB)=[];
end
end
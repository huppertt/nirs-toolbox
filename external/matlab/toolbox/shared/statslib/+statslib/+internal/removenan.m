function [badin,wasnan,varargout]=removenan(varargin)
%REMOVENAN Remove NaN values from inputs
%   [~,WASNAN,X1] = statslib.internal.removenan(Y1) removes missing values
%   from Y1 and returns it as X1. WASNAN is a logical column vector
%   indicating where there were NaN values removed. Y1 is a column vector
%   or a matrix. The type of Y1 can be:
%
%   Categorical       - X1 is categorical, undefined values represents
%                       missing values.
%   Double            - X1 is double. NaN values represents missing values.
%   Single            - X1 is single. NaN values represents missing values.
%   Character matrix  - X1 is a character matrix. Space represents missing
%                       values.
%   Cell              - X1 is a cell array. empty string '' represents
%                       missing values.
%
%  [BAD,WASNAN,X1,x2,...] = statslib.internal.removenan(Y1,Y2,...) accepts any
%  number of input variables Y1,Y2,Y3,..., removes missing values, and
%  returns them as X1, X2,... respectively. All variables should have the
%  same number of observations, and BAD is 0 if that is true. If they to
%  not have the same number, BAD returns the index of a variable whose
%  length does not match the others.
%
%  See also INSERTNAN.

%   Copyright 1993-2014 The MathWorks, Inc.


badin = 0;
wasnan = 0;
n = -1;

% Find NaN, check length, and store outputs temporarily
varargout = cell(nargout,1);
for j=1:nargin
   y = varargin{j};
   if (size(y,1)==1) && (n~=1) 
       y =  y';
   end

   ny = size(y,1);
   if (n==-1)
      n = ny;
   elseif (n~=ny && ny~=0)
      if (badin==0), badin = j; end
   end
   
   varargout{j} = y;

   if (badin==0 && ny>0)
       wasnan = wasnan | any(isnan(y),2);
   end
end

if (badin>0), return; end

% Fix outputs
if (any(wasnan))
   t = ~wasnan;
   for j=1:nargin
      y = varargout{j};
      if (length(y)>0), varargout{j} = y(t,:); end
   end
end

function h = scatterplot(ax,x,y,i,c,m,msize)
%SCATTERPLOT   Scatter plot grouped by index vector.
%   SCATTERPLOT(AX,X,Y,I,C,M,msize) displays a scatter plot of X vs. Y grouped
%   by the index vector I in the axes AX.  
%
%   No error checking.  Use GSCATTER instead.
%
%   Returns handles h to the creates lines. They are arranged so that even if X
%   or Y is a matrix, the first ni elements of hh(:) represent different
%   groups, so they are suitable for use in drawing a legend.
%
%   See also ISCATTER, GSCATTER, GPLOTMATRIX.

%   Copyright 1993-2014 The MathWorks, Inc.

holdState = get(gca,'Nextplot');
holdStateCleanup = onCleanup(@()set(gca,'Nextplot',holdState));

xcols = size(x,2);
ycols = size(y,2);
ncols = max(xcols,ycols);

ni = max(i); % number of groups
if (isempty(ni))
   i = ones(size(x,1),1);
   ni = 1;
end

nm = length(m);
ns = length(msize);
if ischar(c) && isvector(c)
    c = c(:);
end
nc = size(c,1);

% Now draw the plot
for j=1:ni
   ii = (i == j);
   nii = sum(ii);
   if nii==0
       % degenerate case, create lines with NaN data
       for k=1:ncols
           h(j,k) = plot(NaN,NaN, ...
               'Parent', ax, ...
               'LineStyle','none', 'Color', c(1+mod(j-1,nc),:), ...
               'Marker', m(1+mod(j-1,nm)), 'MarkerSize', msize(1+mod(j-1,ns)));
           hold on
       end
   elseif nii==1 && xcols>1 && ycols>1
       % degenerate case, avoid plotting two row vectors in one line
       for k=1:ncols
           h(j,k) = plot(x(ni,k),y(ni,k), ...
               'Parent', ax, ...
               'LineStyle','none', 'Color', c(1+mod(j-1,nc),:), ...
               'Marker', m(1+mod(j-1,nm)), 'MarkerSize', msize(1+mod(j-1,ns)));
           hold on
       end
   else
       % normal case, create lines for each column pair at once       
       h(j,:) = plot(x(ii,:), y(ii,:), ...
           'Parent', ax, ...
           'LineStyle','none', 'Color', c(1+mod(j-1,nc),:), ...
           'Marker', m(1+mod(j-1,nm)), 'MarkerSize', msize(1+mod(j-1,ns))); 
       hold on
    end
end


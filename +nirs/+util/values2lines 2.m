function [colors,lineStyles]=values2lines(v,vrange)

if nargin < 2 || isempty(vrange)
    vmax    = max(abs(v(:)));
    vrange  = vmax*[-1 1];
end

[~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',128) )');
z = linspace(vrange(1), vrange(2), size(cmap,1))';

idx = dsearchn(z,v');

colors = cmap(idx, :);

% line styles
lineStyles = {};
for i = 1:length(idx)
   lineStyles(i,:) = {'LineStyle', '-', 'LineWidth', 4};
end


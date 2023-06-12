function drawon_probe(probe,v,vrange)
% This function will draw the probe using the colors based on the data 

% range to show
if nargin < 3 || isempty(vrange)
    vmax    = max(abs(v(:)));
    vrange  = vmax*[-1 1];
end


% unique types

utypes = unique(probe.link.type, 'stable');


% colormap
[~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',128) )');
z = linspace(vrange(1), vrange(2), size(cmap,1))';

figure;
for iType=1:length(utypes)
    h=subplot(length(utypes),1,iType);
    title(utypes(iType));
    lst = find(ismember(probe.link.type, utypes(iType)));
    vals = v(lst);
    idx = bsxfun(@minus, vals, z);
    [~, idx] = min(abs(idx), [], 1);
    colors = cmap(idx, :);
    
    probe.draw(colors,[],h);
end

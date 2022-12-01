function ChannelSummaryPlot(data)

if(length(data)>1)
    for i=1:length(data)
        nirs.reports.ChannelSummaryPlot(data(i));
    end
    return
end


figure;

ml=data.probe.link;
ml.type=[];
uml=unique(ml);
nchan = size(uml,1)+1;

set(gcf,'Units','normalized');
s=subplot(nchan,1,1);
set(s,'position',[.05 .05+.9-.9*(1)/(nchan+1) .9 1/(nchan+1)]); data.draw([]); legend off;
for i=1:size(uml,1)
    lst=find(ml.source==uml.source(i) & ml.detector==uml.detector(i));
    s=subplot(nchan,1,i+1);
    set(s,'position',[.05 .05+.9-.9*(i+1)/(nchan+1) .9 1/(nchan+1)])
    plot(data.time,data.data(:,lst));
    axis tight;
end


function h=draw_table(tbl,va)
% This function will draw the stats from a table (e.g. as returned by
% ChannelStats.table or nirs.util.roiaverage

if(nargin<2 || isempty(va))
    va={'Beta','beta'};
end


if(ismember('ROI',tbl.Properties.VariableNames))
    rois=tbl.ROI;
    type=rois;
    for i=1:length(type)
        l=strfind(rois{i},':');
        rois{i}(l:end)=[];
        type{i}(1:l)=[];
    end
    tbl.ROI=rois;
    tbl=[tbl table(type)];
    n=length(unique(tbl.Contrast));
    tbl=sortrows(tbl,{'Contrast','type'});
    Names1=unique(tbl.Contrast);
    for i=1:height(tbl)
        Names2{i}=['ROI' tbl.ROI{i} ':' num2str(tbl.fileIdx(i))];
    end
    Names2=unique(Names2);
else
   tbl=sortrows(tbl,{'cond','type'});
   Names1=unique(tbl.cond);
   n=length(unique(tbl.cond));
    for i=1:height(tbl)
        Names2{i}=['Src' num2str(tbl.source(i)) '-Det' num2str(tbl.detector(i))];
    end
    Names2=unique(Names2);
end

types=unique(tbl.type);

X=tbl.(tbl.Properties.VariableNames{find(ismember(tbl.Properties.VariableNames,va))});
if(ismember(va,{'beta','Beta'}))
    E=tbl.(tbl.Properties.VariableNames{find(ismember(tbl.Properties.VariableNames,{'SE','StdErr','se'}))});
else
    E=ones(size(X));
end

figure;
for i=1:length(types)
    lst=ismember(tbl.type,types{i});
    e=E(lst);
    x=X(lst);
    e=reshape(e,[],n);
    x=reshape(x,[],n);
    
    subplot(1,length(types),i);
    nirs.util.bar_err(e,x);
    set(legend(Names1),'Interpreter','none');
    set(title(types{i}),'Interpreter','none');
    set(gca,'XTickLabel',Names2);
     xt=get(gca,'XTick');
    xt=linspace(xt(1),xt(end),length(Names2));
    set(gca,'XTick',xt);
    set(gca,'XtickLabelRotation',45);
end

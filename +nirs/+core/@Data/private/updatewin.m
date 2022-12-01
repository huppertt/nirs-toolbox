function updatewin(varargin)

handles=get(findobj('tag','nirsviewer'),'Userdata');
    l=findobj('type','line','parent',handles.axis_main,'visible','on','Tag','dataline');
    lstim=findobj('type','line','parent',handles.axis_main,'visible','on','Tag','');

rangey=[-1 1];


for idx=1:length(l)
    rangey(idx,1)=min(get(l(idx),'Ydata'));
    rangey(idx,2)=max(get(l(idx),'Ydata'));
    rangex(idx,1)=min(get(l(idx),'Xdata'));
    rangex(idx,2)=max(get(l(idx),'Xdata'));
end

minY=min(rangey(:));
rangeY = max(rangey(:))-minY;

for idx=1:length(lstim)
    yd=get(lstim(idx),'yData');
    yd=yd-min(yd);
    yd=(yd./max(yd))*rangeY*.1+minY-rangeY*.15;
    set(lstim(idx),'yData',yd);
    rangex(length(l)+idx,1)=min(get(lstim(idx),'Xdata'));
   rangex(length(l)+idx,2)=max(get(lstim(idx),'Xdata'));
end
minY=minY-rangeY*.15;

set(handles.axes_mainplot,'Ylim',[minY max(rangey(:))],'Xlim',[min(rangex(:)) max(rangex(:))]);


return


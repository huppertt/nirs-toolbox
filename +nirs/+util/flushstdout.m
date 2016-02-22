function flushstdout

desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
cmdwin = desktop.getClient('Command Window');
cmdwinview = cmdwin.getComponent(0).getViewport.getComponent(0);

while(1)
    
    s=getText(cmdwinview);
    str=s.toCharArray';
    lst=find(double(str)==10);
    
    if(isempty(lst) || double(str(end))==10)
        break
    end
    pause(.1);
    fprintf('\b')
end



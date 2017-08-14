function flushstdout(nlines)

if(nargin==0)
    nlines=12;
end

desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
cmdwin = desktop.getClient('Command Window');
cmdwinview = cmdwin.getComponent(0).getViewport.getComponent(0);

s=getText(cmdwinview);
str=s.toCharArray';
lst=find(double(str)==10);
nlines=min(nlines,length(lst)-1);


for i=1:length(str)-lst(end-nlines)
    fprintf(1,'\b');
end
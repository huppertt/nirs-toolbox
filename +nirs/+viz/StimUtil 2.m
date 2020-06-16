function data=StimUtil(data,showdata)

if(nargin<2)
    showdata=false;
end

%This is a wrapper for the nirs_viewer GUI.
data2=StimUtil_GUI(data,showdata);

%if pressed cancel or closed window then will return empty
if(~isempty(data2))
    data=data2;
end

return

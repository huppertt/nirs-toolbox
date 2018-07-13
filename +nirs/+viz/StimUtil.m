function data=StimUtil(data)
%This is a wrapper for the nirs_viewer GUI.
data2=StimUtil_GUI(data);

%if pressed cancel or closed window then will return empty
if(~isempty(data2))
    data=data2;
end

return

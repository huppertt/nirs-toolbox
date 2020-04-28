<<<<<<< HEAD
function data=StimUtil(data)
%This is a wrapper for the nirs_viewer GUI.
data2=StimUtil_GUI(data);
=======
function data=StimUtil(data,showdata)

if(nargin<2)
    showdata=false;
end

%This is a wrapper for the nirs_viewer GUI.
data2=StimUtil_GUI(data,showdata);
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f

%if pressed cancel or closed window then will return empty
if(~isempty(data2))
    data=data2;
end

return

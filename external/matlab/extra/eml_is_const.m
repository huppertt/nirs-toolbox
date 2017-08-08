function flag=eml_is_const(name)
%dummy function to bypass internal matlab checking 

flag=true;
try
    flag=builtin('eml_is_const',name);
end


return;
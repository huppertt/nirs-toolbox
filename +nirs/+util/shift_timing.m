function data = shift_timing(data,shift)

data.time=data.time+shift;

for i=1:data.stimulus.count
    key=data.stimulus.keys{i};
    st=data.stimulus(key);
    if(isa(st,'nirs.design.StimulusEvents'))
        st.onset=st.onset+shift;
    else
        st.time=st.time+shift;
    end
    data.stimulus(key)=st;
end

if(isa(data,'nirs.core.Data'))
    for i=1:data.auxillary.count
        key=data.auxillary.keys{i};
        aux=data.auxillary(key);
        if(isa(aux,'nirs.core.GenericData'))
            aux.time=aux.time+shift;
            data.auxillary(key)=aux;
        end
    end
end





return
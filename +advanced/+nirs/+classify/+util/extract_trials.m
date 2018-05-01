function varargout=extract_trials(data,type,basis)
% This function extracts single trial estimates
% of the ERP or hemodynamic response

if(nargin<2 || isempty(type))
    type='AR-IRLS';
end
if(nargin<3 || isempty(basis))
    basis=Dictionary;
    
     if(~isempty(strfind(class(data(1).probe),'nirs')))
            basis('default')=nirs.design.basis.Canonical;
     elseif(~isempty(strfind(class(data(1).probe),'eeg')))
            basis('default')=eeg.design.basis.ERP;
     else
         warning('unsupported data type');
         return
     end
end

if(~isa(basis,'Dictionary'))
    b=Dictionary;
    b('default')=basis;
    basis=b;
    clear b;
end

    

Sall=[];

for i=1:length(data)
    
    for sIdx=1:data(i).stimulus.count
        key=data(i).stimulus.keys{sIdx};
        st=data(i).stimulus(key);
        d=data(i);
        for s=1:length(st.onset)
            d(s)=data(i);
            stl=nirs.design.StimulusEvents;
            stl.name=st.name;
            stl.onset=st.onset(s);
            stl.dur=st.dur(s);
            stl.amp=st.amp(s);
            d(s).stimulus=Dictionary;
            d(s).stimulus(key)=stl;
        end
        j=nirs.modules.TrimBaseline;
        j.preBaseline=10;
        j.postBaseline=20;
        d=j.run(d);
        j=nirs.modules.GLM;
        j.type=type;
        j.basis=basis;
        S=j.run(d);
        for j=1:length(S)
            S(j).demographics=d.demographics;
            S(j).demographics('trial')=j;
            S(j).demographics('file')=i;
            S(j).demographics('condition')=key;
        end
        if(isempty(Sall))
            Sall=S(:);
        else
            Sall=[Sall(:); S(:)];
        end
    end
end

varargout{1}=Sall;
if(nargout>1)
    if(~isempty(strfind(class(Sall(1).probe),'nirs')))
        varargout{2}=Sall.HRF;
    elseif(~isempty(strfind(class(Sall(1).probe),'eeg')))
        varargout{2}=Sall.ERP;
    else
        warning('unsupported data type');
        varargout{2}=[];
    end
        
end

return
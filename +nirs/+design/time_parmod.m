function st=time_parmod(stim,time,center,reset)
% This function creates a time-paramtricaly modulated stimulus

if(nargin<3)
    center=true;
end
if(nargin<4)
    reset=false;
end

st=nirs.design.StimulusVector;
st.name=[stim.name '_pmod'];

dt=.5;
t=time';
s=zeros(size(t));

if(reset)
    for i=1:length(stim.onset)
        lst=find(t>=stim.onset(i) & t<=stim.onset(i)+stim.dur(i));
        if(center)
            cnt=lst-mean(lst);
        else
            cnt=lst-lst(1);
        end
        cnt=cnt/max(cnt);
        s(lst)=cnt;
    end
else
    lst=[];
    for i=1:length(stim.onset)
        lst=[lst find(t>=stim.onset(i) & t<=stim.onset(i)+stim.dur(i))];
    end
    if(center)
            cnt=lst-mean(lst);
        else
            cnt=lst-lst(1);
        end
        cnt=cnt/max(cnt);
        s(lst)=cnt;
end

st.time=t;
st.vector=s;
st.convolve_by_default=true;
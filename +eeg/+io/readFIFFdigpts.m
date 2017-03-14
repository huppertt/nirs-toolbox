function digpts = readFIFFdigpts(filename)

[fid,tree]=fiff_open(filename);
[info]=fiff_read_meas_info(fid,tree);
fclose(fid);

% fix this part
T=info.dev_head_t.trans;

id=struct;
for i=1:length(info.dig)
    switch(info.dig(i).kind)
        case(1)
            k='FID';
        case(2)
            k='HPI';
        case(3)
            k='EEG';
        case(4)
            k='extra';
        otherwise
            k=num2str(info.dig(i).kind);
    end
    id.kind{i,1}=k;
    id.X(i,1)=double(info.dig(i).r(1))*1000;
    id.Y(i,1)=double(info.dig(i).r(2))*1000;
    id.Z(i,1)=double(info.dig(i).r(3))*1000;
    id.Frame(i,1)=info.dig(i).coord_frame;
    
    if(info.dig(i).kind==1 & abs(id.X(i))<5)
        id.kind{i,1}='nas';
    elseif(info.dig(i).kind==1 & abs(id.Y(i))<5 & id.X(i)<0)
        id.kind{i,1}='lpa';
    elseif(info.dig(i).kind==1 & abs(id.Y(i))<5 & id.X(i)>0)
        id.kind{i,1}='rpa';
    end
    
end

digpts=unique(struct2table(id),'stable');


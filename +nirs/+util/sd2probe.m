function probe = sd2probe( SD )

if(~isfield(SD,'MeasList') & isfield(SD,'ml'))
    SD.MeasList=SD.ml;
    SD=rmfield(SD,'ml');
end

if(~isfield(SD,'Lambda') & isfield(SD,'lambda'))
    SD.Lambda=SD.lambda;
    SD=rmfield(SD,'lambda');
end

if(~isfield(SD,'SrcPos') & isfield(SD,'srcpos'))
    SD.SrcPos=SD.srcpos;
    SD=rmfield(SD,'srcpos');
end
if(~isfield(SD,'DetPos') & isfield(SD,'detpos'))
    SD.DetPos=SD.detpos;
    SD=rmfield(SD,'detpos');
end
if(~isfield(SD,'AnchorList') & isfield(SD,'al'))
    SD.AnchorList=SD.al;
    SD=rmfield(SD,'al');
end
if(~isfield(SD,'SpringList') & isfield(SD,'sl'))
    SD.SpringList=SD.sl;
    SD=rmfield(SD,'sl');
end
if(~isfield(SD,'DummyPos') & isfield(SD,'dummypos'))
    SD.DummyPos=SD.dummypos;
    SD=rmfield(SD,'dummypos');
end

if(~isfield(SD,'SpatialUnit'))
    SD.SpatialUnit='mm';
end


iSrc    = SD.MeasList(:,1);
iDet    = SD.MeasList(:,2);
wl      = SD.MeasList(:,4);

if ~isfield(SD,'Lambda')
    SD.Lambda = [690 830];
    warning('Assuming wavelengths = [690 830]')
end


wl = SD.Lambda( wl ); wl = wl(:);

link  = table(iSrc,iDet, wl,'VariableNames',{'source','detector','type'});

if(isfield(SD,'SpatialUnit') && strcmp(SD.SpatialUnit,'mm'))
    sc=1;
    
else
    sc=10;  %assume given in cm
    SD.SpatialUnit='cm';
end

probe = nirs.core.Probe( SD.SrcPos*sc, SD.DetPos*sc, link );

%If anchor infomation is present, add it
if(isfield(SD,'AnchorList') && ~isempty(SD.AnchorList))
    PosAll=[SD.SrcPos; SD.DetPos; SD.DummyPos];
    
    cnt=1;
    
    for idx=1:size(SD.AnchorList,1)
        
        
        if(any((SD.SpringList(:,1)==SD.AnchorList{idx,1} | ...
                SD.SpringList(:,2)==SD.AnchorList{idx,1}) & SD.SpringList(:,3)>0))
            Type{cnt,1}='FID-anchor';
            Names{cnt,1}=SD.AnchorList{idx,2};
            
            Units{cnt,1}=SD.SpatialUnit;
            Pos(cnt,:)=PosAll(SD.AnchorList{idx,1},:);
            cnt=cnt+1;
        end
        if(any((SD.SpringList(:,1)==SD.AnchorList{idx,1} | ...
                SD.SpringList(:,2)==SD.AnchorList{idx,1}) & SD.SpringList(:,3)<0))
            Type{cnt,1}='FID-attractor';
            Names{cnt,1}=SD.AnchorList{idx,2};
            
            Units{cnt,1}=SD.SpatialUnit;
            Pos(cnt,:)=PosAll(SD.AnchorList{idx,1},:);
            cnt=cnt+1;
        end
        
    end
    
    
    tbl=table(Names,Pos(:,1),Pos(:,2),Pos(:,3),Type,Units,...
        'VariableNames',probe.optodes.Properties.VariableNames);
    
    probe.optodes=[probe.optodes; tbl];
    
    
end



end


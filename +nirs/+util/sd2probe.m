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

if(isfield(SD,'optpos_reg') && ~isempty(SD.optpos_reg))
    SD.SrcPos=SD.optpos(1:size(SD.SrcPos,1),:);
    SD.DetPos=SD.optpos(1+size(SD.SrcPos,1):end,:);
    SD.SrcPos(:,3)=0;
    SD.DetPos(:,3)=0;
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

% If the distances are all less the 20, then the probe was probably in cm
% and it is mislabeled.  Fix to mm
if(median(probe.distances)<20)
    SD.SpatialUnit='cm';
    sc=10;  %assume given in cm
    probe = nirs.core.Probe( SD.SrcPos*sc, SD.DetPos*sc, link );
end

%If anchor infomation is present, add it
if(isfield(SD,'AnchorList') && ~isempty(SD.AnchorList))
    PosAll=[SD.SrcPos; SD.DetPos; SD.DummyPos];
    Names={};
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
    
    if(~isempty(Names))
        tbl=table(Names,Pos(:,1),Pos(:,2),Pos(:,3),Type,Units,...
            'VariableNames',probe.optodes.Properties.VariableNames);
        
        probe.optodes=[probe.optodes; tbl];
    end
    
    if(isfield(SD,'optpos_reg') && ~isempty(SD.optpos_reg))
            probe1020=nirs.core.Probe1020;
            probe1020.link=probe.link;
            probe1020.optodes=probe.optodes;
            
            T=eye(3);
            if(strcmp(SD.orientation(1),'L'))
                T(1,1)=1;
            end
            if(strcmp(SD.orientation(2),'I'))
                T(2,2)=-1;
            end
            if(strcmp(SD.orientation(3),'P'))
                T(3,3)=-1;
            end
            SP=SD.optpos_reg(1:size(SD.SrcPos,1),:);
            DP=SD.optpos_reg(size(SD.SrcPos,1)+1:end,:);
            SP=SP-ones(size(SP,1),1)*SD.center;
            DP=DP-ones(size(DP,1),1)*SD.center;
            SP=SP*T;
            DP=DP*T;
            
            p2=nirs.core.Probe( SP(:,[1 3 2]), DP(:,[1 3 2]), link );
            
            
            probe1020.optodes_registered=p2.optodes;
            probe=probe1020;
            
    
    end
    
    
end



end


function sd = probe2sd( probe )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        
    assert( isnumeric(probe.link.type) )
    
    [~,~,iWL] = unique( probe.link.type );
    
    sd.nSrcs = size( probe.srcPos, 1);
    sd.nDets = size( probe.detPos, 1);
    sd.Description = ['This is a ' num2str(sd.nSrcs) 'x' num2str(sd.nDets) ' probe.'];
    
    sd.MeasList = ones( length(iWL), 5 );
    sd.MeasList(:,1) = probe.link.source;
    sd.MeasList(:,2) = probe.link.detector;
    sd.MeasList(:,4) = iWL;
    
    sd.MeasListAct=ones(length(iWL),1);
    
    sd.SrcPos = probe.srcPos;
    sd.DetPos = probe.detPos;
    sd.Lambda = unique(probe.link.type,'stable')';
    
end


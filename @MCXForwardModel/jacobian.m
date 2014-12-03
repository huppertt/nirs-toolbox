function [J,meas] = jacobian( obj )
%JACOBIAN Summary of this function goes here
%   Detailed explanation goes here

    origProbe = obj.probe;
    obj.probe = obj.probe.makeUniqueProbe();

    if ~obj.probe.isValid 
        error( 'Probe is not valid.' )
    end

    obj.saveFluence();
    
    meas = load([obj.directory filesep 'measurement.mat']);
    meas = meas.meas;
    
    J.mua = complex( zeros(size(obj.probe.link,1),obj.nLayers,'single') );
    J.kappa = complex( zeros(size(obj.probe.link,1),obj.nLayers,'single') );
    for iLink = 1:size(obj.probe.link,1)
        iDet = obj.probe.link(iLink,2);
        iSrc = obj.probe.link(iLink,1);
        
        src = load([obj.directory filesep 'src' num2str(iSrc) '.mat'],'fluence');
        det = load([obj.directory filesep 'det' num2str(iDet) '.mat'],'fluence');
        
        Jmua = src.fluence .* det.fluence;
        
        [Sx, Sy, Sz] = gradient(src.fluence);
        [Dx, Dy, Dz] = gradient(det.fluence);

%         Sx = filter([1/2 0  -1/2],1,src.fluence,[],1);
%         Sy = filter([1/2 0  -1/2],1,src.fluence,[],2);
%         Sz = filter([1/2 0  -1/2],1,src.fluence,[],3);
% 
%         Dx = filter([1/2 0  -1/2],1,det.fluence,[],1);
%         Dy = filter([1/2 0  -1/2],1,det.fluence,[],2);
%         Dz = filter([1/2 0  -1/2],1,det.fluence,[],3);
        
        Jkappa = ( Sx(:).*Dx(:) + Sy(:).*Dy(:) + Sz(:).*Dz(:) );
        
        for iLayer = 1:obj.nLayers
            mask = obj.image.volume == iLayer;
            J.mua(iLink,iLayer) = - sum( Jmua(mask) ) / meas.data(iLink);
            J.kappa(iLink,iLayer) = - sum( Jkappa(mask) ) / meas.data(iLink);
        end
    end
    
    if obj.cleanup == 1
        rmdir(obj.directory,'s');
    end
end



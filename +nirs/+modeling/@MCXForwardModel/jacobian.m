function [J,meas] = jacobian( obj )
%JACOBIAN Summary of this function goes here
%   Detailed explanation goes here

    % TODO change fft to properly due a dtft
    % for now, choose proper timeStep
    fbin = obj.Fm*1e6 * obj.timeStep * obj.nTimeGates;
    assert( fbin - fix(fbin) < 1e-9 )
    fbin = fix(fbin);
    
    % get a probe in proper format
    obj.probe = obj.probe.makeUniqueProbe();

    assert( obj.probe.isValid() )

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
        
        Jkappa = ( Sx(:).*Dx(:) + Sy(:).*Dy(:) + Sz(:).*Dz(:) );
save('/home/barker/background_figures/jacobian/mcx_jacobian.mat','Jmua','Jkappa','meas');
        for iLayer = 1:obj.nLayers
            mask = obj.image.vol == iLayer;
            J.mua(iLink,iLayer) = - sum( Jmua(mask) ) / meas.data(iLink);
            J.kappa(iLink,iLayer) = - sum( Jkappa(mask) ) / meas.data(iLink);
        end
    end
    
    if obj.cleanup == 1
        rmdir(obj.directory,'s');
    end
end



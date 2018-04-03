function meas = timeResolvedMeas( obj )
%MEASUREMENT Summary of this function goes here
%   Detailed explanation goes here

    % TODO change fft to properly due a dtft
    % for now, choose proper timeStep
%     fbin = obj.Fm*1e6 * obj.timeStep * obj.nTimeGates;
%     assert( fbin - fix(fbin) < 1e-9 )
%     fbin = fix(fbin);
    
    % get a probe in proper format
    probe = obj.probe.makeUniqueProbe();
    
    % either simulate all dets once or all sources once
    % swap if their are less dets than srcs
    if size( probe.detPos,1 ) < size( probe.srcPos,1 )
        probe = probe.swapSD();
    end
    
    % copy fwdModel object with proper probe
    thisObj = obj;
    thisObj.probe = probe;
    
    % check that probe is valid
    assert( probe.isValidSrcs() )
    
    % preallocation
    data = complex( zeros(1,size(probe.link,1),'single') );
    
    % for each source compute flux
    for iSrc = 1:size(probe.srcPos,1)
        cfg = thisObj.getConfig( iSrc,'source' );
        
        tic,
        [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
        for iRep = 2:obj.nRepetitions
            [~,tmp,~] = evalc('mcxlab(cfg)');
            flux.data = (iRep-1)/iRep * flux.data + tmp.data/iRep;
        end
        t = toc;
        disp(['Completed forward model for source ' num2str(iSrc) ' of ' num2str(size(probe.srcPos,1)) ' in ' num2str(t) ' seconds.'])

        flux.data = flux.data;
%         flux.data = fft( flux.data,[],4 );
%         
%         fluence = flux.data(:,:,:,fbin+1);
% if obj.gpuId == 1
% save('/home/barker/neonates/machine_learning/debug.mat')
% end
        lst = find( probe.link(:,1) == iSrc );
        
        pos = probe.detPos( probe.link(lst,2),:);
        pos = pos / thisObj.image.dim(1) + repmat(thisObj.image.origin,[size(pos,1) 1]);
        
        for i = 1:size(flux.data,4)
            tmp = nirs.utilities.median_interp(pos,flux.data(:,:,:,i),obj.image.vol,1);
            tmp( isnan(tmp) ) = 0;
            data(i,lst) = tmp/thisObj.nPhotons/obj.image.dim(1)^3;
        end
            
%         data(lst) = exp(...
%             nirs.utilities.quad_interp(pos,log(fluence),3)...
%             );
        
    end
    
    meas = nirs.Data( data, obj.probe, 1/obj.timeStep, 0, 'MCX' );
    
end


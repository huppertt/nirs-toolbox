function saveFluence( obj )
    
    if ~exist( obj.directory,'dir' )
        mkdir(obj.directory);
    end

 	origProbe = obj.probe;
    newProbe = obj.probe.makeUniqueProbe();
    obj.probe = newProbe;
    
    save( [obj.directory filesep 'probe.mat'],'newProbe','origProbe' );

    if ~obj.probe.isValid 
        error( 'Probe is not valid.' )
    end
    
    Fs = 1/obj.timeStep;
    
	data1 = complex( zeros(1,size(obj.probe.link,1),'single') );
    data2 = complex( zeros(1,size(obj.probe.link,1),'single') );
    
    for iSrc = 1:size(obj.probe.srcPos,1)
        cfg = obj.getConfig( iSrc,'source' );
        
%         tic,
%         [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
%         t = toc;
        
        tic,
        [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
        for iRep = 2:obj.nRepetitions
            [~,tmp,~] = evalc('mcxlab(cfg)');
            flux.data = (iRep-1)/iRep * flux.data + tmp.data/iRep;
        end
        t = toc;
% save( [obj.directory filesep 'src' num2str(iSrc) '_flux.mat'], 'flux' );
        disp(['Completed forward model for source ' num2str(iSrc) ' of ' num2str(size(obj.probe.srcPos,1)) ' in ' num2str(t) ' seconds.'])

        flux.data = flux.data/obj.nPhotons/obj.image.dim(1)^3;%/obj.nRepetitions;
        flux.data = fft( flux.data,[],4 );
        
        fluence = flux.data(:,:,:,fix(obj.Fm*1e6/(Fs/2/obj.nTimeGates))+1);

        save( [obj.directory filesep 'src' num2str(iSrc) '.mat'], 'fluence' );
        
        lst = find( obj.probe.link(:,1) == iSrc );
        
        pos = obj.probe.detPos( obj.probe.link(lst,2),: );
        pos = pos / obj.image.dim(1) + repmat(obj.image.origin,[size(pos,1) 1]);

%         thisR = nirs2.utilities.quad_interp(pos,real(log(fluence)),3);
%         thisI = nirs2.utilities.quad_interp(pos,imag(log(fluence)),3);
        data1(lst) = exp(...
            nirs2.utilities.quad_interp(pos,log(fluence),3)...
            );

    end
    
    for iDet = 1:size(obj.probe.detPos,1)
        cfg = obj.getConfig( iDet,'detector' );
        
%         tic,
%         [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
%         t = toc;
        
        tic,
        [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
        for iRep = 2:obj.nRepetitions
            [~,tmp,~] = evalc('mcxlab(cfg)');
            flux.data = (iRep-1)/iRep * flux.data + tmp.data/iRep;
        end
        t = toc;
        
        disp(['Completed forward model for detector ' num2str(iDet) ' of ' num2str(size(obj.probe.detPos,1)) ' in ' num2str(t) ' seconds.'])

        flux.data = flux.data/obj.nPhotons/obj.image.dim(1)^3;%/obj.nRepetitions;
        flux.data = fft( flux.data,[],4 );
        
        fluence = flux.data(:,:,:,fix(obj.Fm*1e6/(Fs/obj.nTimeGates))+1);

        save( [obj.directory filesep 'det' num2str(iDet) '.mat'], 'fluence' );

        lst = find( obj.probe.link(:,2) == iDet );
        
        pos = obj.probe.srcPos( obj.probe.link(lst,1),: );
        pos = pos / obj.image.dim(1) + repmat(obj.image.origin,[size(pos,1) 1]);

%         thisR = nirs2.utilities.quad_interp(pos,real(log(fluence)),3);
%         thisI = nirs2.utilities.quad_interp(pos,imag(log(fluence)),3);
        data2(lst) = exp(...
            nirs2.utilities.quad_interp(pos,log(fluence),3)...
            );
    end
    
    data = data1/2 + data2/2;
    
    meas = nirs2.Data( data,origProbe,0,obj.Fm,'MCX Forward Model Simulation');
    
	save( [obj.directory filesep 'measurement.mat'], 'meas' );

end


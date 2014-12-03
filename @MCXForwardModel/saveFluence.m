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
        
        tic,
        [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
        t = toc;
        
        disp(['Completed forward model for source ' num2str(iSrc) ' of ' num2str(size(obj.probe.srcPos,1)) ' in ' num2str(t) ' seconds.'])

        flux.data = flux.data/obj.nPhotons;%/obj.nRepetitions;
        flux.data = fft( flux.data,[],4 );
        
        fluence = flux.data(:,:,:,fix(obj.modFreq/(Fs/obj.nTimeGates))+1);

        save( [obj.directory filesep 'src' num2str(iSrc) '.mat'], 'fluence' );
        
        lst = find( obj.probe.link(:,1) == iSrc );
        pos = obj.probe.detPos( obj.probe.link(lst,2),: );
%         data1(lst) = fluence(  sub2ind( size(fluence),...
%             fix( pos(:,1) ),...
%             fix( pos(:,2) ),...
%             fix( pos(:,3) ) ));
        data1(lst) = obj.getMeasFromFluence(pos,fluence);

    end
    
    for iDet = 1:size(obj.probe.detPos,1)
        cfg = obj.getConfig( iDet,'detector' );
        
        tic,
        [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
        t = toc;
        
        disp(['Completed forward model for detector ' num2str(iDet) ' of ' num2str(size(obj.probe.detPos,1)) ' in ' num2str(t) ' seconds.'])

        flux.data = flux.data/obj.nPhotons;%/obj.nRepetitions;
        flux.data = fft( flux.data,[],4 );
        
        fluence = flux.data(:,:,:,fix(obj.modFreq/(Fs/obj.nTimeGates))+1);

        save( [obj.directory filesep 'det' num2str(iDet) '.mat'], 'fluence' );

        lst = find( obj.probe.link(:,2) == iDet );
        pos = obj.probe.srcPos( obj.probe.link(lst,1),: );
%         data2(lst) = fluence(  sub2ind( size(fluence),...
%             fix( pos(:,1) ),...
%             fix( pos(:,2) ),...
%             fix( pos(:,3) ) ));
        data2(lst) = obj.getMeasFromFluence(pos,fluence);
    end
    
    data = data1/2 + data2/2;
%     data = data1;
    
%     if size(newProbe.detPos,1) > size(newProbe.srcPos,1)
%         data = data1;
%     else
%         data = data2;
%     end
    
    meas = nirs.Data( 0,data,origProbe,obj.modFreq,'MCX Forward Model Simulation');
    
	save( [obj.directory filesep 'measurement.mat'], 'meas' );

end


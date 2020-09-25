function saveFluence( obj )
    
    if ~exist( obj.directory,'dir' )
        mkdir(obj.directory);
    end

 	origProbe = obj.probe;
    newProbe = obj.probe; %.makeUniqueProbe();
    obj.probe = newProbe;
    
    save( [obj.directory filesep 'probe.mat'],'newProbe','origProbe' );

%     if ~obj.probe.isValid 
%         error( 'Probe is not valid.' )
%     end
%     
    Fs = 1/obj.timeStep;
    
	data1 = complex( zeros(1,size(obj.probe.link,1),'single') );
    data2 = complex( zeros(1,size(obj.probe.link,1),'single') );
    
    types=obj.probe.types;
    for iW=1:length(types)
        for iSrc = 1:size(obj.probe.srcPos,1)
            cfg = obj.getConfig( [iSrc types(iW)],'source' );
            
            %         tic,
            %         [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
            %         t = toc;
            
            tic,
         
            cfg.nphoton=100;
            obj.nRepetitions=100;
            [~,flux,~,cfg] = evalc('mmclab(cfg)'); % stupid way to suppress mcxlab output
            cfg.nphoton=5E5;
            for iRep = 2:obj.nRepetitions
                disp(['Iter ' num2str(iRep)]);
                tmp=mmclab(cfg);
                %[~,tmp,~] = evalc('mmclab(cfg)');
                flux.data = (iRep-1)/iRep * flux.data + tmp.data/iRep;
            end
            t = toc;
            % save( [obj.directory filesep 'src' num2str(iSrc) '_flux.mat'], 'flux' );
            disp(['Completed forward model for source ' num2str(iSrc) ' of ' num2str(size(obj.probe.srcPos,1)) ' in ' num2str(t) ' seconds.'])
            
            flux.data = flux.data/obj.nPhotons;%/obj.nRepetitions;
            flux.data = fft( flux.data,[],2 );
            
            fluence = flux.data(:,fix(obj.Fm*1e6/(Fs/2/obj.nTimeGates))+1);
            
            save( [obj.directory filesep 'src' num2str(iSrc) '_' num2str(types(iW)) 'nm.mat'], 'fluence' );
            
            lst = find( obj.probe.link.source == iSrc & obj.probe.link.type==types(iW));
            
            pos = obj.probe.detPos( obj.probe.link.detector(lst),: );
            
            %         thisR = nirs.utilities.quad_interp(pos,real(log(fluence)),3);
            %         thisI = nirs.utilities.quad_interp(pos,imag(log(fluence)),3);
            k=dsearchn(cfg.node,pos);
            data1(lst) = fluence(k);
            
        end
        
        for iDet = 1:size(obj.probe.detPos,1)
            cfg = obj.getConfig( [iDet types(iW)],'detector' );
           
            tic,
            [~,flux,~] = evalc('mmclab(cfg)'); % stupid way to suppress mcxlab output
            for iRep = 2:obj.nRepetitions
                [~,tmp,~] = evalc('mmclab(cfg)');
                flux.data = (iRep-1)/iRep * flux.data + tmp.data/iRep;
            end
            t = toc;
            
            disp(['Completed forward model for detector ' num2str(iDet) ' of ' num2str(size(obj.probe.detPos,1)) ' in ' num2str(t) ' seconds.'])
            
            flux.data = flux.data/obj.nPhotons;%/obj.nRepetitions;
            flux.data = fft( flux.data,[],2 );
            
            fluence = flux.data(:,fix(obj.Fm*1e6/(Fs/obj.nTimeGates))+1);
            
            save( [obj.directory filesep 'det' num2str(iDet) '_' num2str(types(iW)) 'nm.mat'], 'fluence' );
            
            lst = find( obj.probe.link.detector == iDet & obj.probe.link.type==types(iW));
            
            pos = obj.probe.srcPos( obj.probe.link.source(lst),: );
            k=dsearchn(cfg.node,pos);
            data2(lst) = fluence(k);
        end
        
        data = data1/2 + data2/2;
    end
    
    meas = nirs.core.Data( data,0,origProbe,obj.Fm);
    
	save( [obj.directory filesep 'measurement.mat'], 'meas' );

end


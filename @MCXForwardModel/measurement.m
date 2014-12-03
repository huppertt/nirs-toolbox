function meas = measurement( obj )
%MEASUREMENT Summary of this function goes here
%   Detailed explanation goes here

    origProbe = obj.probe;
    obj.probe = obj.probe.makeUniqueProbe();
    
    if size( obj.probe.detPos,1 ) < size( obj.probe.srcPos,1 )
        obj.probe = obj.probe.swapSD();
    end
    
    if ~obj.probe.isValidSrcs()
        error( 'Probe is not valid.' )
    end
    
    Fs = 1/obj.timeStep;
    
    data = complex( zeros(1,size(obj.probe.link,1),'single') );
    for iSrc = 1:size(obj.probe.srcPos,1)
        cfg = obj.getConfig( iSrc,'source' );
        
        tic,
        [~,flux,~] = evalc('mcxlab(cfg)'); % stupid way to suppress mcxlab output
        t = toc;
        
        disp(['Completed forward model for source ' num2str(iSrc) ' of ' num2str(size(obj.probe.srcPos,1)) ' in ' num2str(t) ' seconds.'])

        flux.data = flux.data/obj.nPhotons;%/obj.nRepetitions;
        flux.data = fft( flux.data,[],4 );
        
        fluence = flux.data(:,:,:,fix(obj.modFreq/(Fs/obj.nTimeGates))+1);
        
        lst = find( obj.probe.link(:,1) == iSrc );
        pos = obj.probe.detPos( obj.probe.link(lst,2),:);
        
%         data(lst) = fluence( sub2ind( size(fluence),...
%             fix( pos(:,1) ),...
%             fix( pos(:,2) ),...
%             fix( pos(:,3) ) ));

%         for iPos = 1:size(pos,1)
%             x = pos(iPos,1) / obj.image.dim(1);
%             y = pos(iPos,2) / obj.image.dim(1);
%             z = pos(iPos,3) / obj.image.dim(1);
%             
%             [X, Y, Z] = meshgrid( ...
%                 max(floor(x-1),1):min(ceil(x+1),size(fluence,1)),...
%                 max(floor(y-1),1):min(ceil(y+1),size(fluence,2)),...
%                 max(floor(z-1),1):min(ceil(z+1),size(fluence,3)) );
%            	
%             V = fluence( sub2ind( size(fluence), X, Y, Z ) );
%             if length(V) == 1
%                 data(lst(iPos)) = V;
%             else
%                 data(lst(iPos)) = interp3( X, Y, Z, V, x, y, z );
%             end
%         end
        data(lst) = obj.getMeasFromFluence(pos,fluence);

    end
    
    meas = nirs.Data( 0,data,origProbe,obj.modFreq,'MCX Forward Model Simulation');
    
end


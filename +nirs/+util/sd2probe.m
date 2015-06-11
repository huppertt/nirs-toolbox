function probe = sd2probe( SD )

        iSrc    = SD.MeasList(:,1);
        iDet    = SD.MeasList(:,2);
        wl      = SD.MeasList(:,4);
        
        if ~isfield(SD,'Lambda')
            SD.Lambda = [690 830];
            warning('Assuming wavelengths = [690 830]')
        end
        
        
        wl = SD.Lambda( wl ); wl = wl(:);
        
        link  = table(iSrc,iDet, wl,'VariableNames',{'source','detector','type'});
        
        probe = nirs.core.Probe( SD.SrcPos*10, SD.DetPos*10, link );

end


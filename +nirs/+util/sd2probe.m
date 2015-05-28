function probe = sd2probe( SD )

        iSrc    = SD.MeasList(:,1);
        iDet    = SD.MeasList(:,2);
        wl      = SD.MeasList(:,4);
        
        wl( wl==1 ) = 690; % ASSUMED
        wl( wl==2 ) = 830; % ASSUMED
        
        link  = table(iSrc,iDet, wl,'VariableNames',{'source','detector','type'});
        
        probe = nirs.core.Probe( SD.SrcPos*10, SD.DetPos*10, link );

end


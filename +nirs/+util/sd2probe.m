function probe = sd2probe( SD )

        iSrc = SD.MeasList(:,1);
        iDet = SD.MeasList(:,2);
        wavelen  = SD.MeasList(:,4);
        
        wavelen( wavelen==1 ) = 690; % ASSUMED
        wavelen( wavelen==2 ) = 830; % ASSUMED
        
        link = table(iSrc,iDet, wavelen,'VariableNames',{'source','detector','type'});
        
        probe = nirs.Probe( SD.SrcPos*10, SD.DetPos*10, link );


end


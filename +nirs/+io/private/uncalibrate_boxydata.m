function data = uncalibrate_boxydata(data);

%data_calibrated = (data_raw+term)*factor

if(data.CALIBRATIONinfo.Waveform_Calibration_Values_APPLIED)
    WF = data.CalibrationValues.WF;
    ntps=size(data.Data.time,1);
    
    Factor = reshape(WF.Factor,[],3);
    Factor(:,3)=1;
    Term = reshape(WF.Term,[],3);
    
    Term(:,3)=Term(:,3)*2*pi/360;
    
    %assume the data has not been reshuffled
    data.Data.AC=data.Data.AC./(Factor(:,1)*ones(1,ntps))-Term(:,1)*ones(1,ntps);
    data.Data.DC=data.Data.DC./(Factor(:,2)*ones(1,ntps))-Term(:,2)*ones(1,ntps);
    data.Data.Phase=data.Data.Phase./(Factor(:,3)*ones(1,ntps))-Term(:,3)*ones(1,ntps);
    
    data.CALIBRATIONinfo.Waveform_Calibration_Values_APPLIED=false;
end

return
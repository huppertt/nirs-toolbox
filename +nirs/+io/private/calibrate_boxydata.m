function data = calibrate_boxydata(data)

if(data.CALIBRATIONinfo.Waveform_Calibration_Values_APPLIED)
    warning('calibration alread applied... removing');
    data = uncalibrate_boxydata(data);
end

WF = data.CalibrationValues.WF;
ntps=size(data.Data.time,1);

Factor = reshape(WF.Factor,[],3);
Factor(:,3)=1;
Term = reshape(WF.Term,[],3);

Term(:,3)=Term(:,3)*2*pi/360;

%assume the data has not been reshuffled
data.Data.AC=(Factor(:,1)*ones(1,ntps)).*(data.Data.AC+Term(:,1)*ones(1,ntps));
data.Data.DC=(Factor(:,2)*ones(1,ntps)).*(data.Data.DC+Term(:,2)*ones(1,ntps));
data.Data.Phase=data.Data.Phase+Term(:,3)*ones(1,ntps);

data.CALIBRATIONinfo.Waveform_Calibration_Values_APPLIED=true;


return
function data=add_missing_BIDS_info(data)

% % Info that BIDS wants (but isn't always in the file)
Defaults.Manufacturer='n/a';
Defaults.ManufacturersModelName='n/a';
Defaults.SoftwareVersions='n/a';
Defaults.DeviceSerialNumber='n/a';
%Defaults.HardwareFilters='n/a';

Defaults.RecordingDuration=@(data)max(data.time);
Defaults.SamplingFrequency=@(data)data.Fs;

Defaults.NIRSChannelCount=@(data)height(data.probe.link);
Defaults.NIRSSourceOptodeCount=@(data)size(data.probe.srcPos,1);
Defaults.NIRSDetectorOptodeCount=@(data)size(data.probe.detPos,1);

Defaults.SourceType='n/a';
Defaults.DetectorType='n/a';
Defaults.ShortChannelCount =@(data)length(find(data.probe.distances<15));

Defaults.InstitutionName='n/a';
Defaults.InstitutionAddress='n/a';
Defaults.InstitutionalDepartmentName='n/a';
Defaults.CapManufacturer='n/a';
Defaults.CapManufacturersModelName='n/a';

%Defaults.HeadCircumference='n/a';
Defaults.SubjectArtefactDescription='n/a';
Defaults.TaskDescription='n/a';
Defaults.Instructions='n/a';
Defaults.CogPOID='n/a';
Defaults.NIRSPlacementScheme='n/a';

flds=fields(Defaults);

for i=1:length(data)
    for fI=1:length(flds);
        if(~data(i).demographics.iskey(flds{fI}))
            if(isa(Defaults.(flds{fI}),'function_handle'))
                data(i).demographics(flds{fI})=feval(Defaults.(flds{fI}),data(i));
            else
                data(i).demographics(flds{fI})=Defaults.(flds{fI});
            end

        end
    end
    try
        short_channel=(data(i).probe.distances<15);
        short_channel=strrep(strrep(cellstr(num2str(short_channel)),'0','false'),'1','true');
        units=repmat({'nm'},height(data(i).probe.link),1);
        name=strrep(cellstr(horzcat(repmat('S',height(data(i).probe.link),1),...
            num2str(data(i).probe.link.source),...
            repmat(':D',height(data(i).probe.link),1),num2str(data(i).probe.link.detector),...
            repmat('_',height(data(i).probe.link),1),num2str(data(i).probe.link.type),...
            repmat('nm',height(data(i).probe.link),1))),' ','');
        wavelength_nominal=data(i).probe.link.type;
        data(i).probe.link=[data(i).probe.link table(short_channel,name,wavelength_nominal,units)];
    end

end
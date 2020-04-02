function WF = calibrate_from_phantom(calibration_file,PhantomProperties)
%This function will calculate the calibration values (see ISS
%documentation) from one (or more) data files collected on a phantom of
%known properties
%
% Inputs:
%   calibration_file = file name or data file
%                      may be cell array of multiples files/data
%  PhantomProperties = Struct of phantom properties
%        PhantomProperties.Lambda - list of wavelengths (any order)
%        PhantomProperties.MUA -  absorption value (one per wavelength)
%        PhantomProperties.MUS -  scattering value (one per wavelength)
%                       must be cell array if multiple files are used

% Initialize data struct to hold model
AC=[];
DC=[];
Phase=[];
Distances=[];

MUA=[];
MUS=[];

if(~iscell(calibration_file))
    calibration_file={calibration_file};
    PhantomProperties={PhantomProperties};
end

% load and concatinate multiple files

for fileIdx=1:length(calibration_file)
    
    if(isstruct(calibration_file{fileIdx}))
        % this is actually the ISS data
        data=calibration_file{fileIdx};
    else
        data=ReadBoxyData(calibration_file{fileIdx});
    end
    
    %remove any calibrations (if applied)
    data = uncalibrate_boxydata(data);
    
    nmeas = size(data.Data.AC,1);
    ntps=size(data.Data.time,1);
    
    %Term is calculated by dark noise (not from a phantom), so keep the
    %term part of the calibration
    WF = data.CalibrationValues.WF;
    Term = reshape(WF.Term,[],3);    
    local_AC=data.Data.AC+Term(:,1)*ones(1,ntps);
    local_DC=data.Data.DC+Term(:,2)*ones(1,ntps);
    local_Phase=data.Data.Phase;
    
    local_mua = zeros(nmeas,ntps);
    local_mus = zeros(nmeas,ntps);
    
    for idx=1:nmeas
        wIdx=find(PhantomProperties{fileIdx}.Lambda==data.SD.Lambda(data.SD.MeasurementLst(idx,4)));
        local_mua(idx,:)=PhantomProperties{fileIdx}.MUA(wIdx);
        local_mus(idx,:)=PhantomProperties{fileIdx}.MUS(wIdx);
        
    end
    
    local_Dist=reshape(data.Distances',1,[]);
    local_Dist=local_Dist'*ones(1,ntps);
    
    AC=horzcat(AC, local_AC);
    DC=horzcat(DC, local_DC);
    Phase=horzcat(Phase, local_Phase);
    Distances=horzcat(Distances, local_Dist);
    
    MUA = horzcat(MUA,local_mua);
    MUS = horzcat(MUS,local_mus);
end

%Now estimate the calibration from the data
FreqFactor=data.SD.ModFreq*1E6*2*pi/2.26100e+10; 

D = 1./(3*MUA + 3*MUS);  %Diffusion coefficient


SlopeDC=-sqrt(MUA./D);
SlopeAC = -(MUA.^2./D.^2+FreqFactor.^2./D.^2).^(.25).*cos(.5*atan(FreqFactor./MUA));
SlopePhase = (MUA.^2./D.^2+FreqFactor.^2./D.^2).^(.25).*sin(.5*atan(FreqFactor./MUA));

% These are the functions for recovering MUA/MUS.  
%MuAFromACDC(SlopeAC,SlopeDC,FreqFactor);% : OK
%MuAFromACPhase(SlopeAC,SlopePhase,FreqFactor);% : OK
%MuAFromDCPhase(SlopeDC,SlopePhase,FreqFactor);% : BAD  [Don't use this]
%MuSFromACPhase(SlopeAC,SlopePhase,FreqFactor); %: OK

%Expected measurements
Udc = exp(SlopeDC.*Distances)./Distances.^2;
Uac = exp(SlopeAC.*Distances)./Distances.^2;
Uph = SlopePhase.*Distances;

% compute the factor for AC and DC
Factor_DC = median(Udc./DC,2);
Factor_AC = median(Uac./AC,2);
Term_Phase=Uph-Phase;

%Now, store the data to the form used by BOXY
Factor=reshape(WF.Factor,[],3);
Factor(:,1)=Factor_AC;
Factor(:,2)=Factor_DC;
Factor(:,3)=1;  %Factor_Phase
Factor=reshape(Factor,[],1)';

Term=reshape(WF.Term,[],3);
Term(:,3)=median(Term_Phase,2)*360/2/pi;
Term=reshape(Term,[],1)';

WF.Factor = Factor;
WF.Term=Term;


return
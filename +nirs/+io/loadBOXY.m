function data = loadBOXY(filenames,CalibrationFile,CalibrationPhantom)
% Inputs:
%   calibration_file = file name or data file
%                      may be cell array of multiples files/data
%  PhantomProperties = Struct of phantom properties
%        PhantomProperties.Lambda - list of wavelengths (any order)
%        PhantomProperties.MUA -  absorption value (one per wavelength)
%        PhantomProperties.MUS -  scattering value (one per wavelength)
%                       must be cell array if multiple files are used
       
warning('nirs.io.loadBOXY.m function has hardcoded wavelength info specific for my ISS system; please edit');

wavelengthlst{1}=[1 16 17 32];  %676
wavelengthlst{2}=[2 15 18 31];  %690
wavelengthlst{3}=[3 14 19 30]; %750
wavelengthlst{4}=[4 13 20 29 8 9 24 25]; %788
wavelengthlst{5}=[5 12 21 28]; %800
wavelengthlst{6}=[6 11 22 27]; %808
wavelengthlst{7}=[7 10 23 26]; %830

lambda=[676 690 750 788 800 808 830];

    if ischar( filenames )
        filenames = {filenames};
    end
    
    if(nargin==3)
        WF = calibrate_from_phantom(CalibrationFile,CalibrationPhantom)
    else
        if(nargin==2)
            warning('Both calibration file (BOXY) and Properties need to be assigned');
        end
        
        WF=[];
    end
    
    data = nirs.core.Data.empty;
    
    % iterate through cell array
    for iFile = 1:length(filenames)
        try
            % load data as a struct
            ISSdata = ReadBoxyData(filenames{iFile},WF,wavelengthlst,lambda);
            
            ml=ISSdata.Data.MeasurementList;
            ml(:,5)=ISSdata.Distances(sub2ind(size(ISSdata.Distances),ml(:,2),ml(:,1)));
            ml(:,6)=ml(:,5);
            ml(:,5)=round(ml(:,5)/.5)*.5;
            [a]=unique(ml(:,[2 5]),'rows');
            ISSdata.SD.NumSrc=size(a,1);
            ISSdata.SD.SrcPos = [10*ones(ISSdata.SD.NumSrc,1) linspace(-40,40,ISSdata.SD.NumSrc)' zeros(ISSdata.SD.NumSrc,1)];
            ISSdata.SD.DetPos = [-10*ones(ISSdata.SD.NumDet,1) linspace(-40,40,ISSdata.SD.NumDet)' zeros(ISSdata.SD.NumDet,1)];
            
            for i=1:size(a,1)
                lst=find(ml(:,2)==a(i,1) & ml(:,5)==a(i,2));
                ml(lst,1)=i;
            end
            ml(:,3)=1;
            ISSdata.SD.MeasList=ml(:,1:4);
                        
            %convert ISSdata to the new class structure
            data(iFile)=nirs.core.Data;
            data(iFile).description=filenames{iFile};
            data(iFile).time=ISSdata.Data.time;
            data(iFile).probe=nirs.util.sd2probe(ISSdata.SD);
            data(iFile).probe.fixeddistances=ml(:,6);
            data(iFile).Fm=ISSdata.SD.ModFreq;
            
%             ac = abs( data.data(:,lst) );
%             phs = angle( data.data(:,lst) );
            ISSdata.Data.Phase=phase_unwrap(ISSdata.Data.Phase,ml(:,6),data(iFile).Fm);
            d=ISSdata.Data.AC.*(cos(ISSdata.Data.Phase)+1i*sin(ISSdata.Data.Phase));
            data(iFile).data=d;
            
          
        catch err
            if(~exist('ISSdata') || ~isfield(ISSdata,'Data') || isempty(ISSdata.Data))
                 disp('Empty file found (skipping):');
                 disp(filenames{iFile});
            else
                warning(err.message)
            end
            
        end
        
        
    end

return



function jobs = write_paper(folder)
% This runs a complete analysis including image reconstruction and saving
% all the files into tiff files.


if(nargin<1)
    folder=fullfile(pwd,['analysis ' datestr(now)]);
else
    if(exist(folder)~=7)
        mkdir(folder);
    end
end

folder2=fullfile(folder,'probe');
folder1=fullfile(folder,'image');
disp(['data will be saved in: ' folder ]);

if(exist(folder2)~=7)
    mkdir(folder2);
end

jobs=nirs.modules.ImportData();
jobs.Input='raw';



jobs=nirs.modules.RunMatlabCode(jobs);
jobs.FunctionHandle=@()evalin('base','close all;');



jobs=nirs.modules.RemoveStimless(jobs);
jobs = nirs.modules.Resample(jobs);
jobs.Fs = 5; % resample to 5 Hz
jobs = nirs.modules.OpticalDensity( jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='dOD';

if(evalin('base','isa(raw(1).probe,''nirs.core.Probe1020'')'))
    disp('adding image recon code');
    if(exist(folder1)~=7)
        mkdir(folder1);
    end
    jobs=nirs.modules.Assert(jobs);
    jobs.throwerror=true;
    jobs.condition=@(data)isa(data.probe,'nirs.core.Probe1020');
    jobs = nirs.modules.TrimBaseline( jobs );
    jobs.preBaseline   = 30;
    jobs.postBaseline  = 30;
    jobs = nirs.modules.AR_IRLS(jobs );
    jobs = nirs.modules.ExportData(jobs);
    jobs.Output='SubjStats_dOD';
    jobs = nirs.modules.ImageReconMFX(jobs);
    
    Slab = nirs.forward.ApproxSlab;
    probe=evalin('base','raw(1).probe');
    Slab.probe=probe.swap_reg;
    lambda=unique(probe.link.type);
    Slab.prop=nirs.media.tissues.brain(.7,50,lambda);
    Slab.mesh=probe.getmesh();
    Slab.mesh=Slab.mesh(end);
    Jac=Slab.jacobian('spectral');
    jobs.probe('default')=Slab.probe;
    jobs.jacobian('default')=Jac;
    jobs.formula='beta ~ -1 + cond';
    jobs.mesh=Slab.mesh;
    jobs.basis=nirs.inverse.basis.gaussian(Slab.mesh,30);
    
    jobs = nirs.modules.ExportData(jobs);
    jobs.Output='ImageStats';
    
    jobs=nirs.modules.RunMatlabCode(jobs);
    jobs.FunctionHandle=@()evalin('base',[...
        'ImageStats.printAll(''tstat'',[],[],[],[],''' folder1 ''',''tiff'')']);
    
    
    
end

jobs=nirs.modules.ImportData(jobs);
jobs.Input='dOD';
jobs.override=true;
jobs = nirs.modules.BeerLambertLaw( jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='Hb';
jobs = nirs.modules.TrimBaseline( jobs );
jobs.preBaseline   = 30;
jobs.postBaseline  = 30;
jobs = nirs.modules.AR_IRLS(jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='SubjStats_Hb';
jobs = nirs.modules.MixedEffects(jobs );
jobs.formula       = 'beta ~ -1 + cond';  % See help fitlme for examples
jobs = nirs.modules.ExportData(jobs);
jobs.Output='GroupStats_Hb';

if(evalin('base','isa(raw(1).probe,''nirs.core.Probe1020'')'))
    
    jobs=nirs.modules.RunMatlabCode(jobs);
    jobs.FunctionHandle=@()evalin('base','ROIs=nirs.util.convertlabels2roi(GroupStats_Hb.probe);');
    jobs=nirs.modules.RunMatlabCode(jobs);
    jobs.FunctionHandle=@()evalin('base','ROItable=nirs.util.roiAverage(GroupStats_Hb,ROIs);');
    
    jobs=nirs.modules.RunMatlabCode(jobs);
    jobs.FunctionHandle=@()evalin('base','disp(ROItable(ROItable.q<=0.05,:));');
    
    jobs=nirs.modules.RunMatlabCode(jobs);
    jobs.FunctionHandle=@()evalin('base','GroupStats_Hb.probe.defaultdrawfcn=''10-20'';');
    
end

jobs=nirs.modules.RunMatlabCode(jobs);
jobs.FunctionHandle=@()evalin('base',[...
    ' GroupStats_Hb.printAll(''tstat'',[],''q<0.05'',''' folder2 ''',''tiff'')']);

jobs=nirs.modules.RunMatlabCode(jobs);
jobs.FunctionHandle=@()evalin('base',['writetable(ROItable,''' folder '' filesep 'ROIStats.txt'')']);

jobs=nirs.modules.RunMatlabCode(jobs);
jobs.FunctionHandle=@()evalin('base',['disp([''Congratulations.  Your paper has been submitted! '' char(9786) '' Just kidding'']);'...
    'disp('' But when you write it, please cite:'')']);


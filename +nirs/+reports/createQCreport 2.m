function tbl=createQCreport(raw)

% make summary of info
tbl=nirs.createDemographicsTable(raw);
d.Subjects = unique(tbl.StudyID);
for i=1:length(d.Subjects)
    lst=find(ismember(tbl.StudyID,d.Subjects{i}));
    d.NumScans(i,1)=length(lst);
    dates={};
    for j=1:length(lst)
        a=dir(raw(lst(j)).description);
        dates={dates{:} datestr(a(1).datenum,'mmm-dd-yyyy')};
    end
    ud=unique(dates);
    for j=1:length(ud)
   % d.Date{i,j}=ud{j};
    end
    for j=1:length(lst)
     d.ScanTimes{i,j}=raw(lst(j)).time(end);      
    end
    o=[]; sni=[]; mm=[];
    for j=1:length(lst)
       stim=nirs.getStimNames(raw(lst(j)));
       ons=0;
       o(j)=0;
       for k=1:length(stim)
            s=raw(lst(j)).stimulus(stim{k});
            ons=ons+length(s.onset);
            o(j)=o(j)+length(s.onset);
        end
         d.Conditions{i,j}=[num2str(ons) ' Events'];
           
    end
   
    lst(find(vertcat(d.ScanTimes{i,:})<60 & o'==0))=[];
     for j=1:length(lst)
         [s,m]=nirs.math.structnoiseindex(raw(lst(j)).data);
         sni=[sni s];
         mm=[mm m];
     end
      d.SNI(i,1)=nanmean(sni);
      d.MotionReject(i,1)=nanmean(mm);
    d.View{i,1}=['<a href="matlab: nirs.viz.QCReport(raw([' num2str(lst') ']))">view</a>'];
end 

tbl=struct2table(d);
sni=[];
function write_phoebe_probe(probe1020,filename)

filename=strtok(filename,'.');  % remove the extension

SD = nirs.util.probe2sd(probe1020);
save([filename '.SD'],'SD','-MAT');


mesh=probe1020.getmesh;
fid=mesh(1).fiducials;

mainpts = fid(ismember(fid.Name,{'nas','rpa','lpa'}),:);
mainpts.Name{ismember(mainpts.Name,'nas')}='nz';
mainpts.Name{ismember(mainpts.Name,'rpa')}='ar';
mainpts.Name{ismember(mainpts.Name,'lpa')}='al';


fid=fopen([filename '.txt'],'w');
for i=1:3
    fprintf(fid,'%s:\t%d\t%d\t%d\n',mainpts.Name{i},mainpts.X(i),mainpts.Y(i),mainpts.Z(i));
end

probe1020=probe1020.swap_reg;
for i=1:size(probe1020.srcPos,1)
     fprintf(fid,'s%d:\t%d\t%d\t%d\n',i,probe1020.srcPos(i,1),probe1020.srcPos(i,2),probe1020.srcPos(i,3));
end
for i=1:size(probe1020.detPos,1)
     fprintf(fid,'d%d:\t%d\t%d\t%d\n',i,probe1020.detPos(i,1),probe1020.detPos(i,2),probe1020.detPos(i,3));
end

fclose(fid);



function probe2xml(probe,filename)

fid = fopen(filename,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" ?>\n');
fprintf(fid,'<system>\n');
fprintf(fid,'<probe>\n');

fprintf(fid,'\t<sources>\n');
for i=1:size(probe.srcPos,1);
    fprintf(fid,'\t\t<srcpos> <index> %d </index><x> %g </x> <y> %g </y> </srcpos>\n',...
        round(i-1),probe.srcPos(i,1),probe.srcPos(i,2));
end

if(isa(probe,'nirs.core.Probe1020'))
    p2=probe.swap_reg;
    
    for i=1:size(p2.srcPos,1);
        fprintf(fid,'\t\t<srcpos3D> <index> %d </index> <x> %g </x> <y> %g </y> <z> %g </z> </srcpos3D>\n',...
            round(i-1),p2.srcPos(i,1),p2.srcPos(i,2),p2.srcPos(i,3));
    end
end
fprintf(fid,'\t</sources>\n\n');


fprintf(fid,'\t<detectors>\n');
for i=1:size(probe.detPos,1);
    fprintf(fid,'\t\t<detpos> <index> %d </index> <x> %g </x> <y> %g </y> </detpos>\n',...
        round(i-1),probe.detPos(i,1),probe.detPos(i,2));
end

if(isa(probe,'nirs.core.Probe1020'))
    p2=probe.swap_reg;
    
    for i=1:size(p2.detPos,1);
        fprintf(fid,'\t\t<detpos3D> <index> %d </index> <x> %g </x> <y> %g </y> <z> %g </z> </detpos3D>\n',...
            round(i-1),p2.detPos(i,1),p2.detPos(i,2),p2.detPos(i,3));
    end
end
fprintf(fid,'\t</detectors>\n\n');

fprintf(fid,'\t<measurements>\n');
ml=round(unique([probe.link.source probe.link.detector],'rows'));
for i=1:size(ml,1)
      fprintf(fid,'\t\t<ml> <src> %d </src> <det> %d </det> </ml>\n',ml(i,1)-1,ml(i,2)-1);
end
fprintf(fid,'\t</measurements>\n');

fprintf(fid,'</probe>');
fprintf(fid,'</system>\n');
fclose(fid);
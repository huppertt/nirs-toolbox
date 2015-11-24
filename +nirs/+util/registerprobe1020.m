function probe1020 = registerprobe1020(probe,headsize)
% This function registers a probe to the 10-20 system

if(nargin>1)
    probe1020 = nirs.core.Probe1020([],headsize);
else
    probe1020 = nirs.core.Probe1020;
end
probe1020.link=probe.link;
probe1020.optodes=probe.optodes;

mesh=probe1020.getmesh;
probe = nirs.util.registerProbe2Mesh(mesh(1),probe);
probe1020.optodes_registered=probe.optodes;





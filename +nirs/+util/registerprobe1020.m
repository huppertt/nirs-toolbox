function probe1020 = registerprobe1020(probe)
% This function registers a probe to the 10-20 system

probe1020 = nirs.core.Probe1020;
probe1020.link=probe.link;
probe1020.optodes=probe.optodes;

mesh=probe1020.getmesh;
probe = nirs.util.registerProbe2Mesh(mesh(1),probe);
probe1020.optodes_registered=probe.optodes;





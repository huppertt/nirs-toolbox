function probe1020 = atlasview2probe1020(file,sym)

if(nargin<2)
    sym=false;
end

% This function will convert an atlasviewer probe to a probe1020
MESH=nirs.forward.NirfastBEM();
MESH=MESH.loadBEM_atlasviewer(file);
MESH.mesh(1).fiducials.Draw(:)=false;
tmp=load(file,'probe');
probe=nirs.util.sd2probe(tmp.probe);

probe1020 = nirs.core.Probe1020;
probe1020.optodes=probe.optodes;
probe1020.optodes_registered=MESH.probe.optodes;
probe1020.link=MESH.probe.link;

MESH.mesh(1).transparency=.1;
MESH.mesh(2).transparency=1;

probe1020=probe1020.regsister_mesh2probe(MESH.mesh);

% Use the corrdinates in the atlasiewer file instead of our re-registration
fidtbl =probe1020.getmesh.fiducials;
fidtbl.Draw=[];
fidtbl(find(ismember(fidtbl.Type,'10-20')),:)=[];
probe1020.optodes_registered=fidtbl;

if(sym)
    probe1020=nirs.registration.realign_symetric(probe1020);
end

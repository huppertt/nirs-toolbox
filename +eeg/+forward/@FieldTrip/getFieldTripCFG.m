function cfg = getFieldTripMeshes( obj )
    probe = obj.probe;
    
    cnt=1;
    for i=1:length(obj.prop)
        if(~isnan(obj.prop(i)))
            bnd(cnt).tri=obj.mesh.mesh(i).faces;
            bnd(cnt).pnt=obj.mesh.mesh(i).nodes; 
            conduct(cnt)=obj.prop(i);
            cnt=cnt+1;
            
        end
    end


    vol.bnd=bnd;
    vol.type='bemcp';
    vol.unit='mm';
    vol.cfg.method='bem';
    vol.cfg.trackconfig='off';
    vol.cfg.checkconfig='loose';
    vol.cfg.checksize=100000;
    vol.cfg.showcallinfo='yes';
    vol.mat=1;
    
    % Conductivities for model.
    vol.cond =conduct;
    
    cnt=1;
    sens=[];
    sens.type='eeg';
    sens.pnt=[];  %For now, just use a subset of the surface as EEG positions;  TODO- pay more attention to avoid the face, neck etc.
    sens.ori=[];
    
    % Register the probe to the head
    [~,lst]=ismember(lower(obj.probe.electrodes.Name),lower(obj.mesh.mesh(1).fiducials.Name));
    if(length(lst)==height(obj.probe.electrodes))
        pts=[obj.mesh.mesh(1).fiducials.X(lst) obj.mesh.mesh(1).fiducials.Y(lst) ...
            obj.mesh.mesh(1).fiducials.Z(lst)];
         [k,d]=dsearchn(bnd(1).pnt,pts);
    end
    for i=1:height(obj.probe.electrodes)
        sens.pnt(i,:)=pts(i,:);
        sens.ori(i,:)=sens.pnt(i,:)/norm(sens.pnt(i,:));
        sens.label{i,1}=obj.probe.electrodes.Name{i};
        sens.chantype{i,1}='eeg';
    end
    
    cfg = [];
    cfg.elec=sens;
    cfg.channel = {'EEG'};   % the used channels; but subtract bad channels, ex: '-MLP31', '-MLO12'
    cfg.grid.pos = obj.mesh.mesh(end).nodes;              % use the White matter as the source points
    cfg.grid.inside = 1:size(obj.mesh.mesh(end).nodes,1); % all source points are inside of the brain
    cfg.grid.unit='mm';
    cfg.vol = vol;  % volume conduction model
    
end


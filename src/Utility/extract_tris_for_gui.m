function extract_tris_for_gui(scnpth, scnpth_raw, scnfls, scnfls_raw)
for n = 1:length(scnfls)
    %---------------------------------------------------------------------------
    % Load the raw data
    [colobj_raw,tri] = read_Heges_ply(scnpth_raw,[scnfls_raw{n}],1);        
    %---------------------------------------------------------------------------
    % Load the processed data
    dat     = load([scnpth,'/',scnfls{n}]);
    %---------------------------------------------------------------------------
    % Crop the raw data and the corresponding triangulation
    colobj = dat.colobj;
    %---------------------------------------------------------------------------
    % Recreate cropping on the raw scan from the
    ris    = find( (sum(double(colobj_raw.Color),2) < 2) | ...
        (colobj_raw.Location(:,1) < colobj.XLimits(1)) | (colobj_raw.Location(:,1) > colobj.XLimits(2)) | ...
        (colobj_raw.Location(:,2) < colobj.YLimits(1)) | (colobj_raw.Location(:,2) > colobj.YLimits(2)) | ...
        (colobj_raw.Location(:,3) < colobj.ZLimits(1)) | (colobj_raw.Location(:,3) > colobj.ZLimits(2)));
    locs        = colobj_raw.Location;
    cols        = colobj_raw.Color;
    [p,t]       = update_tris_due_to_rmvnodes(locs,tri,ris);
    
    name=scnfls{n}(1:end-4);
    eval(['save ',scnpth,'/',name,'_rawtris.mat p t'])
    size(t)
    colobj
    %----------------------------------------------------------------------
    figure
    trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','cyan','linestyle','none')
    camlight left
    drawnow
    pause(0.1)
end
end



%-------------------------------------------------------------------------------
% 
% Load all the iphone scanning data
% 
%-------------------------------------------------------------------------------
function allscans = load_iphone_scan_dat_and_rawtris_hr(scnpth,scnfls,scnpth_raw,scnfls_raw,dbg_flg)

%-------------------------------------------------------------------------------
for n = 1:length(scnfls)
    %---------------------------------------------------------------------------
    % Load the raw data
    [colobj_raw,tri_raw,locs_raw,cols_raw] = read_Heges_ply(scnpth_raw,[scnfls_raw{n}],0);        
    %---------------------------------------------------------------------------
    % Load the processed data
    dat     = load([scnpth,'/',scnfls{n}]);
    %---------------------------------------------------------------------------
    % Crop the raw data and the corresponding triangulation
    colobj = dat.colobj;
    %     %---------------------------------------------------------------------------
    %     % Recreate cropping on the raw scan from the
    %     ris    = find( (sum(double(colobj_raw.Color),2) < 2) | ...
    %         (colobj_raw.Location(:,1) < colobj.XLimits(1)) | (colobj_raw.Location(:,1) > colobj.XLimits(2)) | ...
    %         (colobj_raw.Location(:,2) < colobj.YLimits(1)) | (colobj_raw.Location(:,2) > colobj.YLimits(2)) | ...
    %         (colobj_raw.Location(:,3) < colobj.ZLimits(1)) | (colobj_raw.Location(:,3) > colobj.ZLimits(2)));
    %     locs        = colobj_raw.Location;
    %     cols        = colobj_raw.Color;
    %     [p,t]       = update_tris_due_to_rmvnodes(locs,tri,ris);
    %     cols(ris,:) = [];
    %     colobj_raw2 = pointCloud(p, 'Color', uint8(cols(:,1:3)));
    %
    %     if dbg_flg == 1
    %         figure
    %         pcshow(colobj_raw2)
    %         figure
    %         trisurf(t(:,1:3),p(:,1),p(:,2),p(:,3),'facecolor','cyan','linestyle','none')
    %     end
    %----------------------------------------------------------
    % Make the raw triangles and electrode found scan data
    % consistent
    [tri,colobj_elf] = make_tri_elfound_consistent(colobj_raw,tri_raw,colobj,0);
    cols_elf         = colobj_elf.Color;
    locs_elf         = colobj_elf.Location;

    %     %----------------------------------------------------------------------
    %     elunlab = dat.newtdats;
    %     el = 90*pi/180;% str2double(get(handles.ini_rx,'String'))*pi/180;
    %     Rx          = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
    %     tdat        = (Rx*(dat.tdat'))';
    %     elunlab     = (Rx*(elunlab'))';
    %     locs        = (Rx*(dat.colobj.Location'))';
    %     colobj      = pointCloud(locs, 'Color', dat.colobj.Color);
    %----------------------------------------------------------------------
    elunlab = dat.newtdats;
    el = 90*pi/180;% str2double(get(handles.ini_rx,'String'))*pi/180;
    Rx          = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
    tdat        = (Rx*(dat.tdat'))';
    elunlab     = (Rx*(elunlab'))';
    locs        = (Rx*(colobj_elf.Location'))';
    colobj_elf  = pointCloud(locs, 'Color', colobj_elf.Color);

    %----------------------------------------------------------------------
    allscans(n).colobj  = colobj_elf;
    allscans(n).tdat    = tdat;
    allscans(n).elunlab = elunlab;
    allscans(n).tri     = tri;
end
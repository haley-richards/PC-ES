%--------------------------------------------------------------------------
%
% Make the raw triangles and electrode found scan data consistent
%
%--------------------------------------------------------------------------
function [t,colobj] = make_tri_elfound_consistent(colobj_raw,tri,colobj,dbg_flg)


%--------------------------------------------------------------------------
% Make sure we only are considering unique nodes in the raw scan
locs_raw                = round(double(colobj_raw.Location),6);
cols_raw                = colobj_raw.Color;
[tri,locs_raw,used_nds] = remove_repeated_nodes(tri,locs_raw);
cols_raw                = cols_raw(used_nds,:);

%--------------------------------------------------------------------------
% Similarly, make sure we are only considering unique nodes in the 
% electrode found scan
locs_elf       = round(double(colobj.Location),6);
cols_elf       = colobj.Color;
[locs_elf,ia0] = unique(locs_elf,'rows');
cols_elf       = cols_elf(ia0,:);

%--------------------------------------------------------------------------
% Keep only the associated nodes between the two sets
[locs_elf,ia2,ib2]  = intersect(locs_elf,locs_raw,'rows');
cols_elf            = cols_elf(ia2,:);
ris                 = setdiff(1:size(locs_raw,1),ib2);
[p,t]               = update_tris_due_to_rmvnodes(locs_raw,tri,ris);

%--------------------------------------------------------------------------
% Return the updated point cloud
colobj       = pointCloud(locs_elf, 'Color', cols_elf);

if dbg_flg == 1
    figure
    trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','cyan','linestyle','none')
    camlight left
    view(3)

    figure
    trisurf(t,locs_elf(:,1),locs_elf(:,2),locs_elf(:,3),'facecolor','cyan','linestyle','none')
    camlight left
    view(3)
end

end
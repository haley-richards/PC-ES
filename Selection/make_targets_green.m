function [pc_out]=make_targets_green(pc_in, tdat, ncol)
pc_out=pc_in;
is_targ=find(isnan(tdat(:, 1))==0);
pts=tdat(is_targ, :);
for i=1:size(pts, 1);
%------------------------------------------------------------------------------------
%convert 2d ui point selections onto 3d point cloud
dist=sqrt((pc_in.Location(:, 1)-pts(i, 1)).^2+(pc_in.Location(:, 2)-pts(i, 2)).^2+(pc_in.Location(:, 3)-pts(i, 3)).^2);
[sorted, idx]=sort(dist);
pc_out.Color(idx(1:ncol), :)=repmat([0, 250, 0], ncol, 1);
end
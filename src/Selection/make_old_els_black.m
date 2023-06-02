function [pc_out]=make_old_els_black(pc_in)
pc_out=pc_in;
cols=pc_in.Color;
is_1= find((cols(:, 1)==250)==1);
is_2=find((cols(:, 2)==0)==1);
is_3=find((cols(:, 3)==0)==1);
i1= intersect(is_1, is_2);    
is_red= intersect(i1, is_3);
for k=1:size(is_red, 1)
cols(is_red(k), 1:3)=[0, 0, 0];
end
pc_out.Color=cols;
% is_targ=find(isnan(tdat(:, 1))==0);
% pts=tdat(is_targ, :);
% for i=1:size(pts, 1);
% %------------------------------------------------------------------------------------
% %convert 2d ui point selections onto 3d point cloud
% dist=sqrt((pc_in.Location(:, 1)-pts(i, 1)).^2+(pc_in.Location(:, 2)-pts(i, 2)).^2+(pc_in.Location(:, 3)-pts(i, 3)).^2);
% [sorted, idx]=sort(dist);
% pc_out.Color(idx(1:ncol), :)=repmat([0, 0, 0], ncol, 1);
end

function [pc_sub_out, proj_data]= project_scan_onto_axes_and_subdivide_v4(colobj, n_part) 
  %colobj=proj_pc_to_xy_plane_v2(colobj);%*do this better
  [model,inliers] = pcfitsphere(colobj,.015);
  locs=colobj.Location;
  cols=colobj.Color;
  colobj_2=pointCloud(locs(inliers, :), 'Color', cols(inliers, :));
  locs_cart=locs(:, :);
  [az, el, r]=cart2sph(locs_cart(:, 1),locs_cart(:, 2),locs_cart(:, 3));
  [idx, C] = kmeans([az, el, r], n_part,'Distance', 'cityblock' );
%   figure;
% plot3(locs_cart(idx==1,1),locs_cart(idx==1,2),locs_cart(idx==1,3), 'r.','MarkerSize',12)
% hold on
% plot3(locs_cart(idx==2,1),locs_cart(idx==2,2),locs_cart(idx==2,3), 'g.','MarkerSize',12)
% plot3(locs_cart(idx==3,1),locs_cart(idx==3,2),locs_cart(idx==3,3), 'b.','MarkerSize',12)
% plot3(locs_cart(idx==4,1),locs_cart(idx==4,2),locs_cart(idx==4,3), 'k.','MarkerSize',12)
% % C_cart=sph2cart(C);
% % plot3(C_cart(:,1),C_cart(:,2),C_cart(:,3),'kx',...
% %      'MarkerSize',15,'LineWidth',3) 
% % legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4', 'Centroids',...
% %        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off
  az_tol=pi/3;
  el_tol=az_tol;
  for i=1:n_part
      in=find(idx==i);
      %get overlap regions
      in_overlap=find(and(abs(az-C(i, 1))<az_tol, abs(el-C(i, 2))<el_tol)==1);
      in_final=unique([in; in_overlap]);
      C_out(i, :)=[mean(az(in_final)), mean(el(in_final)), mean(r(in_final))];
      pc_sub_out{i}=(in_final);
  end

  proj_data= C;
  %proj_data(:, 3)=repmat(mean(r), 4, 1);
 
  
  % ---------------------------------------------------------
  
end

 
 
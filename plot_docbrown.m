function plot_docbrown(mdl,el_col)
hold on
plot3(mdl.sens_ps(:,1),mdl.sens_ps(:,2),mdl.sens_ps(:,3),['.',el_col],'markersize',24)
% plot3(mdl.lmrk_ps(:,1),mdl.lmrk_ps(:,2),mdl.lmrk_ps(:,3),'.r','markersize',24)
%--- plot edges
cnxs = mdl.cnxs;
nds  = mdl.sens_ps; % [mdl.sens_ps; mdl.lmrk_ps];
size(nds)
plot3( ...
    [nds(cnxs(:,1),1) nds(cnxs(:,2),1)]', ... 
    [nds(cnxs(:,1),2) nds(cnxs(:,2),2)]', ...
    [nds(cnxs(:,1),3) nds(cnxs(:,2),3)]', ...
    '-b')
view(3)
axis equal 
axis off
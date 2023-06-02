function [newtdats0, colobj_out] = run_single_prj_alg_v2_sub(newtdats0,handles, colobj0)
% clear
% clc
% close all
% load testdat
dbg_flg = 0;

%-------------------------------------------------------------------------------
[colobj,tdat,eclks,cent0] = update_PLY_pts(handles);
colobj_sub=colobj;
%-------------------------------------------------------------------------------
% Rotate the PLY file and points found to the current frame
if isempty(newtdats0) 
    newtdats = [];
else    
    newtdats = rotate_general_pts(newtdats0,handles);
end
% Get new points
[tdats, colobj_out]   = run_prj_alg_v2_sub(colobj_sub,newtdats,dbg_flg, colobj0);
if size(tdats,1) > 0
    tdats     = convert_pts_back(tdats,handles,cent0);
    % Add the points (in the original frame) to the full set of points
    newtdats0         = [newtdats0; tdats];
end
drawnow

% %-------------------------------------------------------------------------------
% % To check that it worked...
% handles.eclks = (1.02*sampdirs0);
% update_ptcld_plot(hObject,handles)

% %-------------------------------------------------------------------------------
% figure
% set(gcf,'position',[258         213        1025         739])
% hold on
% az2        = str2double(get(handles.cur_az,'String'))*pi/180;
% el2        = str2double(get(handles.cur_el,'String'))*pi/180;
% [colobj,tdat,eclks,cent0] = update_PLY_pts(handles,[az2 el2]);
% newtdats = rotate_general_pts(newtdats0,handles,[az2 el2]);
% MS = 100;
% pcshow(colobj,'Markersize',MS);
% plot3(1.01*newtdats(:,1),1.01*newtdats(:,2),1.01*newtdats(:,3),'.m','markersize',20)
% axis equal
% view(2)

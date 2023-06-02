function [tdat_out, colobj_out] = run_single_targ_alg_v2_sub(tdat0,handles, colobj0)
% clear
% clc
% close all
% load testdat
dbg_flg = 0;
tdat_out=tdat0;
%-------------------------------------------------------------------------------
[colobj,tdat,eclks,cent0] = update_PLY_pts(handles);
colobj_sub=colobj;
%-------------------------------------------------------------------------------
% Rotate the PLY file and points found to the current frame
if isempty(tdat0) 
    tdats = [];
else    
    tdats = rotate_general_pts(tdat0,handles);
end
% Get new points
[newtdats, targ_label, colobj_out]  = run_targ_alg_v2_sub(colobj_sub,tdats, dbg_flg, colobj0);

%tdat_out0(targ_label, :)=newtdats;

if size(tdats,1) > 0
    tdats1    = convert_pts_back(tdats,handles,cent0);
    tdats2    =convert_pts_back(newtdats, handles, cent0);
    % Add the points (in the original frame) to the full set of points
    tdat_out(targ_label, :)=tdats2;
end
drawnow
%tdat_out=tdat;
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
end
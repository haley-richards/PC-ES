%-------------------------------------------------------------------------------
%
% Run the projection algorithm
%
%-------------------------------------------------------------------------------
function [labtdats,tdats0] = run_singleiview_labeling_prj_alg(handles)
% clear
% clc
% close all
% load testdat
dbg_flg = 0;

%-------------------------------------------------------------------------------
tdats0    = handles.tdat;
newtdats0 = handles.newtdats;
if isfield(handles,'labtdats') == 0
    labtdats  = NaN*newtdats0(:,1);
else
    labtdats  = handles.labtdats;
    if size(labtdats,1) < size(newtdats0,1)
        % Need to add blank labels
        disp('Need to add blank labels')
        labtdats = [labtdats; NaN*ones(size(newtdats0,1)-size(labtdats,1),1)];
    end
end
sclf_ps   = [];

%-------------------------------------------------------------------------------
% Load the nominal cap
nomdat = load_nominal_cap(0);

%-------------------------------------------------------------------------------
% Rotate the PLY file and points found to the current frame
[colobj,tdat,eclks,cent0] = update_PLY_pts(handles);

%---------------------------------------------------------------------------
% Rotate the points to the current view
if size(newtdats0,1) > 0
    newtdats = rotate_general_pts(newtdats0,handles);
else
    newtdats = [];
end

%---------------------------------------------------------------------------
% Get the best-fit nominal electrode positions
[Bout,elns,cnxs,sclf_ps] = bestfit_nomel_noplot(tdat,nomdat,colobj,sclf_ps);
% figure
% hold on
% plot3(tdat(:,1),tdat(:,2),tdat(:,3),'.k','markersize',18)
% plot3(Bout(:,1),Bout(:,2),Bout(:,3),'.r','markersize',18)
%---------------------------------------------------------------------------
% Get new points
[labtdats]   = run_prj_alg_labeling(colobj,newtdats,Bout,elns,cnxs,labtdats,dbg_flg);
% if dbg_flg == 1
drawnow

%---------------------------------------------------------------------------
% Update the list of labeled targets
%     invelns = zeros(257,1);
%     for k = 1:length(elns)
%         invelns(elns(k)) = k;
%     end
efounds = find( isnan(labtdats) == 0);
for k = 1:length(efounds)
    tdats0(labtdats(efounds(k)),1:3) = newtdats0(efounds(k),1:3);
end



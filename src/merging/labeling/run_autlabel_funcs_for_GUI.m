%--------------------------------------------------------------------------
%
% Run the projection algorithm
%   1. Perform a rigid transformation of the point cloud to the nominal cap
%
%--------------------------------------------------------------------------
function [tdat,elunlab] = run_autlabel_funcs_for_GUI(handles)
% clear
% clc
% close all
% load testdat
% % save testdat
% % clear
dbg_flg = 0;
addpath(genpath('S:/digihisto/Ethan/EKM_utility'))
handles
%--------------------------------------------------------------------------
% Load the nominal cap
nomdat = load_nominal_cap(0);

% norm( handles.newtdats - handles.eclks)
% error('stop')
%--------------------------------------------------------------------------
[colobj,tdat,eclks,cent0] = update_PLY_pts(handles);
allscans0.colobj  = colobj;
allscans0.tdat    = tdat;
allscans0.eclks   = eclks;
allscans0.elunlab = eclks;
 
%--------------------------------------------------------------------------
[allscans,nomdat] = bestfit_scan_to_nominal(allscans0,nomdat,0);



%--------------------------------------------------------------------------
dethr   = 10/1000;       % Maximum distance to consider for automatic 
                        % labeling of electrodes (m)
tdat    = allscans.tdat;
colobj  = allscans.colobj;
elunlab = allscans.elunlab;
[tdat,n_unlab,n_ellab,elunlab] = auto_label_func(elunlab,tdat,colobj,nomdat,dethr,0);
set(handles.message_texts,'String',[num2str(size(elunlab,1)),' unlabeled, ',num2str(length(find( (isnan(tdat(:,1)) == 0) ))),' labeled electrodes'])


%--------------------------------------------------------------------------
% Convert the unlabeled electrodes back to the original frame (not best
% fitting to the nominal data
T       = allscans.T;
A       = T(1:3,1:3);
t       = T(1:3,4);
tdat    = (A'*(tdat' - repmat(t,1,size(tdat,1))))';
tdat    = convert_pts_back(tdat,handles,cent0);
elunlab = (A'*(elunlab' - repmat(t,1,size(elunlab,1))))';
elunlab = convert_pts_back(elunlab,handles,cent0);





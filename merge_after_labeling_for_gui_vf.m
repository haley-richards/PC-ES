%-------------------------------------------------------------------------------
% 
% merge scans and autolabel
% 
%-------------------------------------------------------------------------------
function [colobj_mrg, tdat_mrg, allscans]= merge_after_labeling_for_gui_vf(scnpth, scnfls, scnpth_raw, scnfls_raw, valpth, valfnm, subj_name, elthick)

% addpath(genpath('S:\digihisto\Ethan\EKM_Utility'))
% addpath(genpath('S:\eit\Haley\stroke\Hedges_Gui\gui_v7'))

elthick   = 12/1000;    % electrode thickness
dethr     = 10/1000;     % Maximum distance to consider for automatic 
                        % labeling of electrodes (m)

 
%-------------------------------------------------------------------------------
% Load the nominal cap
nomdat = load_nominal_cap(0);
%------------------------------------------------------------------------------
% Load the photogrammetry data
pgdat = load_photogram_dat(valpth,valfnm);
%-------------------------------------------------------------------------------
% Load all the iphone scanning data
allscans0 = load_iphone_scan_dat_and_rawtris(scnpth,scnfls,scnpth_raw,scnfls_raw,1);
close all
save([subj_name, '_allscans0.mat'], 'allscans0')
%get rid of unlabeled points that are too close to other points
dist_thresh=.005;
count=1;
err_cnt=0;
for scn=1:length(allscans0)
    count=1;
    elunlab=[];
    elunlab0=allscans0(scn).elunlab;
    n_el_before=size(elunlab0, 1);
    disp(['before= ', num2str(n_el_before)]);
    for j=1:size(elunlab0, 1)
        el_test=elunlab0(j, :);
     dist=vecnorm(elunlab0-el_test, 2, 2);
     dist_sort=sort(dist);
     min_dist=dist_sort(2);
     if min_dist<dist_thresh
      disp('pair too close!');
      err_cnt=err_cnt+1;
     else
         elunlab(count, :)=elunlab0(j, :);
          count=count+1;
     end
    end
     n_el_after=size(elunlab, 1);
    disp(['after= ', num2str(n_el_after)]);
   allscans0(scn).elunlab=elunlab;
   
end

%--------------------------------------------------------------------------
% Loop through all the scans and construct labels
lab_res = zeros(length(allscans0),4);
for n = 1:length(allscans0)
    disp(['Auto label: scan ',num2str(n),' of ',num2str(length(allscans0))])
    n_unlab_ini = size(allscans0(n).elunlab,1);
    n_ellab_ini = length(find( (isnan(allscans0(n).tdat(:,1)) == 0)));
    [tdat,n_unlab,n_ellab] = auto_label_func_hr(allscans0(n).elunlab,allscans0(n).tdat,allscans0(n).colobj,nomdat,dethr,0);
    allscans0(n).tdat = tdat;
    lab_res(n,:) = [n_unlab_ini n_ellab_ini n_unlab n_ellab];
  %allscans0(n).tdat=allscans0(n).tdat(1:256, :);
end
disp(' ')
disp('*******')
disp(['                  Number of scans: ',num2str(length(allscans0))])
disp(['Initial Ave. unlabeled electrodes: ',num2str(round(mean(lab_res(:,1)),1)),'+/-',num2str(round(std(lab_res(:,1)),1))])
disp(['                  Percent Labeled: ',num2str(round(mean(lab_res(:,4)./lab_res(:,1)*100),1)),'+/-',num2str(round(std(lab_res(:,4)./lab_res(:,1)*100),1))])
disp(' ')


%-------------------------------------------------------------------------------
% Simulaneously merge all scans into the nominal cap
[allscans,nomdat] = simulmerg2nominal(allscans0,nomdat,0);

save([subj_name, '_allscans.mat'], 'allscans')
%--------------------------------------------------------------------------
% Compute the normal of each electrode
allscans2 = comp_elec_nvecs(allscans,elthick,0);

% Merge scan by merging labeled electrodes by from most intersections to
% least.dd
[~,tdat_mrg_temp] = merge_scans_for_target_error_finding(allscans2,scnfls, 1)
%arrange allscans from most to least electrodes for better merging
for n_scn=1:length(allscans2)
    n_elunlab(n_scn)=size(allscans2(n_scn).elunlab, 1);
end
[k, idx]=sort(n_elunlab,'descend');
allscans_sort=allscans2;
for n_scn=1:length(allscans2)
    allscans_sort(n_scn)=allscans2(idx(n_scn));
end
[colobj_mrg,tdat_mrg] = merge_labeled_scans_tdatdirect_v2(allscans2,0);
colobj_mrg2 = pointCloud(colobj_mrg.Location, 'Color', colobj_mrg.Color);% ,'Normal',normals);
%--------------------------------------------------------------------------
%infer missing electrodes 
% [tdat_mrg_fullset] = find_missing_elec_labels_err(nomdat,pgdat, tdat_mrg,colobj_mrg,1)
%--------------------------------------------------------------------------
% Calculate the overall error
shrt_name=subj_name;
[valtrans, tdat_mrg_filt, rms_with_filt, max_err, n_el] = calc_valerr_v2(pgdat,tdat_mrg,colobj_mrg,1, subj_name);
[tdat_mrg_fullset, rms_fullset] = find_missing_elec_labels_err(nomdat,pgdat, tdat_mrg_filt,colobj_mrg,0, subj_name);
%nomps0 = calc_valerr_v2(nomdat,tdat_mrg_fullset,colobj_mrg,1,[scnpth_raw,'/',subj_name]);

%--------------------------------------------------------------------------
% Save the merged data
el = -90*pi/180;% str2double(get(handles.ini_rx,'String'))*pi/180;
Rx        = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
locs      = double(colobj_mrg.Location);
tdat_mrg  = (Rx*(tdat_mrg'))';
locs      = (Rx*(locs'))';
colobj_mrg= pointCloud(locs, 'Color', colobj_mrg.Color);
tdat      = tdat_mrg;
colobj    = colobj_mrg;
rms_s     = rms_with_filter;
rms_fullset=rms_fullset
file_ini  = 'merged scan';

eval(['save ',scnpth,'/merged_scan file_ini colobj', ...
        ' tdat'])
end


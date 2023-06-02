%-------------------------------------------------------------------------------
% 
% merge scans and autolabel
% 
%-------------------------------------------------------------------------------
function [colobj_mrg, tdat_mrg, allscans]= merge_after_labeling_for_gui_vf_no_pg(scnpth, scnfls, scnpth_raw, scnfls_raw, subj_name, elthick, savepth)
%no photogrammetry data involved

elthick   = elthick/1000;    % electrode thickness
dethr     = 10/1000;     % Maximum distance to consider for automatic 
                        % labeling of electrodes (m)

%-------------------------------------------------------------------------------
% Load the nominal cap
nomdat = load_nominal_cap(0);
%-------------------------------------------------------------------------------
% Load all the iphone scanning data
allscans0 = load_iphone_scan_dat_and_rawtris_hr(scnpth,scnfls,scnpth_raw,scnfls_raw,1);
%close all
eval(['save ',savepth,'/allscans0.mat allscans0'])
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

eval(['save ',savepth,'/allscans.mat allscans'])
%--------------------------------------------------------------------------
% Compute the normal of each electrode
allscans2 = comp_elec_nvecs(allscans,elthick,0);

% Merge scan by merging labeled electrodes by from most intersections to
% least.dd
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
%save the image
saveas(gcf,[savepth,'/',subj_name,'.png'])
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
newtdats  = [];
file_ini  = 'merged scan';

eval(['save ',savepth,'/merged_scan file_ini colobj', ...
        ' tdat'])
end


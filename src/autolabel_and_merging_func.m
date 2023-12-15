%--------------------------------------------------------------------------
%
% Autolabel the scans and then merge them.
% The steps performed here are the following:
%   1. Get rid of too closely spaced unlabeled electrodes (if ones are too
%      close, we just average them together.
%   2. Autolabeling: Loop through all the scans and find electrode labels
%   3. Transform each iphone scan into the nominal EEG frame
%   4. Compute the normal of each electrode and project the point down to
%      the scalp based on the electode thickness
%   5. Merge scan by merging labeled electrodes by from most intersections to
%      least.
%   6. Calculate errors
%      a. RMS of unmerged electrodes found with an outlier threshold
%      b. RMS of merged electrodes found with an outlier threshold
%      c. RMS of merged electrodes found after filling in missing
%         electrodes
%
%--------------------------------------------------------------------------
function outinf = autolabel_and_merging_func( ...
    allscans0, nomdat, elthick, dethr, dfill, outthresh, dbg_flg, pgdat )


%--------------------------------------------------------------------------
% electrode thickness
elthick=elthick/1000;   % Convert the electrode thickness to m from mm

%--------------------------------------------------------------------------
%   1. Get rid of too closely spaced unlabeled electrodes (if ones are too
%      close, we just average them together.
% Get rid of unlabeled points that are too close to other points. We loop
% through the set of unlabeled electrodes and check if any two points are
% too close together. If they are too close together, then should we
% average them? - let's go for this - but then make sure we have a unique
% set of points at the end of the loop.
dist_thresh = 0.005;
for scn = 1:length(allscans0)
    count       = 1;
    elunlab     = [];
    elunlab0    = allscans0(scn).elunlab;
    n_el_before = size(elunlab0, 1);
    %----------------------------------------------------------------------
    if dbg_flg == 1
        disp(['before= ', num2str(n_el_before)]);
    end
    for j=1:size(elunlab0, 1)
        el_test         = elunlab0(j, :);
        dist            = vecnorm(elunlab0-el_test, 2, 2);
        [dist_sort,sid] = sort(dist);
        iless           = find(dist_sort < dist_thresh);
        if length(iless) == 1 % No points too close
            elunlab(count,:) = elunlab0(j, :);
            count            = count+1;
        else
            % Take the average of these points
            disp('pair too close!');
            elunlab(count,:) = mean(elunlab0(sid(iless), :),1);
            count            = count+1;
        end
    end
    %----------------------------------------------------------------------
    elunlab                = unique(elunlab,'rows');
    n_el_after             = size(elunlab, 1);
    allscans0(scn).elunlab = elunlab;
    if dbg_flg == 1
        disp(['after= ', num2str(n_el_after)]);
    end
end

%--------------------------------------------------------------------------
% 2. Autolabeling: Loop through all the scans and find electrode labels
autolab_start = toc;
lab_res = zeros(length(allscans0),4);
for n = 1:length(allscans0)
    %----------------------------------------------------------------------
    if dbg_flg == 1
        disp(['Auto label: scan ',num2str(n),' of ',num2str(length(allscans0))])
    end

    %----------------------------------------------------------------------
    % Run the autolabel algorithm
    [tdat,n_unlab,n_ellab] = auto_label_func_hr_loc(allscans0(n).elunlab,allscans0(n).tdat,allscans0(n).colobj,nomdat,dethr,0);
    allscans0(n).tdat = tdat;

    %----------------------------------------------------------------------
    % Record the
    lab_res(n,:) = [ ...
        size(allscans0(n).elunlab,1) ...    % Initial number of unlabeled electrodes
        length(find( (~isnan(allscans0(n).tdat(:,1)) ))) ... % Initial number of labeled electrodes
        n_unlab ...                     % Final number of unlabeled electrodes
        n_ellab];                       % Final number of labeled electrodes
    % allscans0(n).tdat=allscans0(n).tdat(1:256, :);
end
if dbg_flg == 1
    disp(' ')
    disp('*******')
    disp(['                  Number of scans: ',num2str(length(allscans0))])
    disp(['Initial Ave. unlabeled electrodes: ',num2str(round(mean(lab_res(:,1)),1)),'+/-',num2str(round(std(lab_res(:,1)),1))])
    disp(['                  Percent Labeled: ',num2str(round(mean(lab_res(:,4)./lab_res(:,1)*100),1)),'+/-',num2str(round(std(lab_res(:,4)./lab_res(:,1)*100),1))])
    disp(' ')
end
autolab_end = toc;

mrge_start = toc;
%--------------------------------------------------------------------------
% 3. Transform each iphone scan into the nominal EEG frame
[allscans,nomdat] = simulmerg2nominal_loc(allscans0,nomdat,dbg_flg);


%--------------------------------------------------------------------------
% 4. Compute the normal of each electrode and project the point down to the
%    scalp based on the electode thickness
allscans2 = comp_elec_nvecs(allscans,elthick,0);

%--------------------------------------------------------------------------
% 5. Merge scan by merging labeled electrodes by from most intersections to
%    least.
[colobj_mrg,tdat_mrg,tri_mrg,nconxs,rmserrs_mrg,omitted_scans] = merge_labeled_scans_tdatdirect_v4_loc(allscans2,0);

mrge_end = toc;

%--------------------------------------------------------------------------
% Get a final best-fit between the labeled electrodes and the nominal cap
[Bout,Tout,bestsclf,rmserrs_mrg(end+1)] = bestfit_nomel_noplot_fullset_loc(tdat_mrg,nomdat,dbg_flg);


%--------------------------------------------------------------------------
tdat_mrg_fullset = find_missing_elec_labels(nomdat, tdat_mrg, colobj_mrg, dfill, dbg_flg);

%--------------------------------------------------------------------------
% Output without errors
outinf.colobj_mrg  = colobj_mrg;
outinf.tdat_mrg    = tdat_mrg;
outinf.tri_mrg     = tri_mrg;
outinf.allscans    = allscans;
outinf.nconxs      = nconxs;
outinf.rmserrs_mrg = rmserrs_mrg;
outinf.time_to_autolab     = autolab_end - autolab_start;    
outinf.time_to_merge       = mrge_end - mrge_start;
outinf.omitted_scans       = omitted_scans;
outinf.tdat_full           = tdat_mrg_fullset;

%--------------------------------------------------------------------------
nargin
if nargin == 8
    %----------------------------------------------------------------------
    % 6. Calculate errors
    %    a. RMS of unmerged electrodes found with an outlier threshold
    %    b. RMS of merged electrodes found with an outlier threshold
    %    c. RMS of merged electrodes found after filling in missing electrodes
    [rms_sep_with_filter, mdn_sep_with_filter,max_sep_with_filter, n_sep_els] = calc_valerr_matlab_pgdat_sep(pgdat,allscans2,outthresh,dbg_flg);
    %--------------------------
    [tdat_mrg_filt, rms_with_filter,mdn, max_err_mrg, n_el,diff_mrg] = calc_valerr_matlab_pgdat_v3(pgdat,tdat_mrg,colobj_mrg,outthresh,dbg_flg);
    %--------------------------
    [tdat_mrg_fullset, rms_fullset,mdn_fullset, max_err_full,diff_full] = find_missing_elec_labels_err_v4( ...
        nomdat, pgdat, tdat_mrg, colobj_mrg, dfill, dbg_flg);


    %----------------------------------------------------------------------
    % Output the errors
    outinf.rms_sep_with_filter = rms_sep_with_filter;
    outinf.mdn_sep_with_filter = mdn_sep_with_filter;
    outinf.max_sep_with_filter = max_sep_with_filter;
    outinf.n_sep_els           = n_sep_els;
    outinf.rms_with_filter     = rms_with_filter;
    outinf.mdn_err             = mdn;
    outinf.max_err_mrg         = max_err_mrg;
    outinf.n_el                = n_el;
    outinf.rms_fullset         = rms_fullset;
    outinf.mdn_fullset         = mdn_fullset;
    outinf.max_err_full        = max_err_full;    
    outinf.diff_mrg            = diff_mrg;
    outinf.diff_full           = diff_full;

end
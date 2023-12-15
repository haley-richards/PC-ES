%--------------------------------------------------------------------------
%
% find_missing_elec_labels_err_v4 - this method uses a distance threshold. 
% The algorithm proceeds as following:
%   1. Find all the missing electrodes
%   2. Calculate the number of electrodes within the distance threshold 
%      for each missing electrode. We do a best-fit using all available 
%      electrodes to get nominal locations of all the electrodes. 
%   3. Iterative fill in missing electrodes, always using the missing 
%      electrode that has the most points within the distance threshold
%      and using a best fit only with points within that threshold
%
%--------------------------------------------------------------------------
function tdat_mrg = find_missing_elec_labels(nomdat, tdat_mrg,colobj_mrg,dfill,dbg_flg)

%--------------------------------------------------------------------------
%   1. Find all the missing electrodes
tdat_mrg=tdat_mrg(1:256, :);
is_fnd = find(~isnan(tdat_mrg(:,1)) ); % Index of found electrodes
is_mis = find( isnan(tdat_mrg(:,1)) ); % Index of missing electrodes

if dbg_flg == 1
    disp([num2str(length(is_fnd)),' Electrodes found'])
end


%--------------------------------------------------------------------------
%   3. Iterative fill in missing electrodes, always using the missing
%      electrode that has the most points within the distance threshold
%      and using a best fit only with points within that threshold
while ~isempty(is_mis)
    %----------------------------------------------------------------------
    %   2. Calculate the number of electrodes within the distance threshold
    %      for each missing electrode. We do a best-fit using all available
    %      electrodes to get nominal locations of all the electrodes.
    % * Bout is the best-fit nominal electrodes on matching points
    % * Cout is the complete nominal set of electrodes
    [Bout,Tout,bestsclf] = bestfit_nomel_noplot_fullset_loc(tdat_mrg,nomdat,0);
    Cout = (Tout*[diag(bestsclf)*(nomdat.eeg_mps'); ones(1,size(nomdat.eeg_mps, 1))])';
    %----------------------------------
    num_ins = zeros(length(is_mis),1);
    for n = 1:length(is_mis)
        is = find( sqrt(sum( (tdat_mrg - repmat(Cout(is_mis(n),1:3),size(tdat_mrg,1),1)).^2,2)) < dfill);
        num_ins(n) = length(is);
    end
    %----------------------------------------------------------------------
    [tmp,i] = max(num_ins);
    is = find( sqrt(sum( (tdat_mrg - repmat(Cout(is_mis(i),1:3),size(tdat_mrg,1),1)).^2,2)) < dfill);
    %----------------------------------------------------------------------
    % Update the best fit, only using points within the distance threshold
    tdat_tmp       = NaN*tdat_mrg;
    tdat_tmp(is,:) = tdat_mrg(is,:);
    [Bout,Tout,bestsclf] = bestfit_nomel_noplot_fullset_loc(tdat_tmp,nomdat,0);
    Cout = (Tout*[diag(bestsclf)*(nomdat.eeg_mps'); ones(1,size(nomdat.eeg_mps, 1))])';

    if dbg_flg == 1
        if rand > 0.5
            figure;hold on
            plot3(tdat_mrg(:,1),tdat_mrg(:,2),tdat_mrg(:,3),'.k','markersize',8)
            plot3(tdat_tmp(:,1),tdat_tmp(:,2),tdat_tmp(:,3),'or','markersize',4)
            plot3(Cout(:,1),Cout(:,2),Cout(:,3),'.b','markersize',8)
            plot3(Cout(is_mis(i),1),Cout(is_mis(i),2),Cout(is_mis(i),3),'.m','markersize',18)
            lbl_fmt_fig('X (m)','Y (m)','','','Z (m)',12)
            legend('Tdat_{merged}','Tdat_{local}','Nominal Local fit','Replacement point')
            clear
        end
    end
    %----------------------------------------------------------------------
    % Fill in the point and remove it from the missing list
    tdat_mrg(is_mis(i),:) = Cout(is_mis(i),1:3);
    is_mis(i) = [];

end


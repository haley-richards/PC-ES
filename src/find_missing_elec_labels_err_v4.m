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
function [valtrans, rms, mdn, max_err, diff] = find_missing_elec_labels_err_v4(nomdat,pgdat, tdat_mrg,colobj_mrg,dfill,dbg_flg)

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

%--------------------------------------------------------------------------
B_full  = [tdat_mrg(1:256, :) ones(256,1)]';
A_full  = [pgdat(1:256, :), ones(256, 1)]';
%-----------------------------------------------------------------------
% Get registered points from each frame
[T, rmserr] = transform_loc(A_full, B_full, 1);
valtrans    = (T*B_full)';
diff        = vecnorm((valtrans(:, 1:3)-pgdat(1:256, :)), 2, 2);

%--------------------------------------------------------------------------
rms     = sqrt(1/length(diff)*sum( diff.^2 ));
mdn     = median(abs(diff(isnan(diff)==0)));
max_err = max(abs(diff));


%--------------------------------------------------------------------------
if dbg_flg == 1
    %     disp(['mean rms= ', num2str(rms*1000), ' mm'])
    %     disp(['median rms= ', num2str(mdn*1000), ' mm'])
    %     disp(['max rms= ', num2str(max_err*1000), ' mm'])
    %     save(['fullset_merging_data_', subj_name], 'mdn', 'rms', 'max_err')
    disp(['tdat full set / photogrammetry RMS Error: ',num2str(rmserrs(k)*1000),' mm'])

    %---------------------------------------------------------------------------
    % Plot the nominal eeg cap and the
    figure
    set(gcf,'position',[680         394        1109         584])
    subplot(1,2,1)
    hold on
    plot3(pgdat(:,1),pgdat(:,2),pgdat(:,3),'.k','markersize',16)
    plot3(valtrans(:,1),valtrans(:,2),valtrans(:,3),'.r','markersize',16)
    for n = 1:size(tdat_mrg,1)
        plot3([tdat_mrg(n,1) valtrans(n,1)],[tdat_mrg(n,2) valtrans(n,2)],[tdat_mrg(n,3) valtrans(n,3)],'-g')
    end
    axis equal
    grid on
    box on
    FS = 12;
    xlabel('X','fontsize',FS,'FontName','times')
    ylabel('Y','fontsize',FS,'FontName','times')
    zlabel('Z','fontsize',FS,'FontName','times')
    set(gca,'FontSize',FS,'FontName','times')
    set(gcf,'color','w');

    axis off

    subplot(1,2,2) % figure
    hold on
    MS = 100;
    pcshow(colobj_mrg,'Markersize',MS);
    axis equal
    grid on
    box on
    FS = 12;
    xlabel('X','fontsize',FS,'FontName','times')
    ylabel('Y','fontsize',FS,'FontName','times')
    zlabel('Z','fontsize',FS,'FontName','times')
    set(gca,'FontSize',FS,'FontName','times')
    set(gcf,'color','w');
    axis off


end
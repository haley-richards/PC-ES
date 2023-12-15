%--------------------------------------------------------------------------
%
% calc_valerr_matlab_pgdat_v2: Calculate quantitative metrics between the
% electrodes found from the iPhone and the photogrammetry results. 
%
%--------------------------------------------------------------------------
function [valtrans, rms_with_filter,mdn, max_err, n_el, diff,rmserr_nom_fit] = calc_valerr_matlab_pgdat_v3(pgdat,tdat_mrg,colobj_mrg,outthresh,dbg_flg)

%--------------------------------------------------------------------------
% Find the indices corresponding to found electrodes
tdat_mrg = tdat_mrg(1:256, :);
is_fnd   = find( ~isnan(tdat_mrg(:,1)) );

if dbg_flg == 1
    disp([num2str(length(is_fnd)),' Electrodes found'])
    figure
    hold on
    plot3(tdat_mrg(is_fnd,1),tdat_mrg(is_fnd,2),tdat_mrg(is_fnd,3),'.k','markersize',16)
    plot3(pgdat(is_fnd,1),pgdat(is_fnd,2),pgdat(is_fnd,3),'.r','markersize',16)
end

%--------------------------------------------------------------------------
% Best-fit the iphone and photogrammetry results. Specifically map the
% iphone to the photogrammetry frame, i.e. T is such that A = T * B.
B = [tdat_mrg(is_fnd,:) ones(length(is_fnd),1)]';
A = [pgdat(is_fnd,:) ones(length(is_fnd),1)]';
[T, rmserr_fit]   = transform_loc(A, B, 1);
if dbg_flg == 1
    disp(['RMS error: ',num2str(rmserr*1000),' mm'])
end

%--------------------------------------------------------------------------
% Calculuate the 2-norm vector between the two scans
valtrans = (T*([tdat_mrg ones(size(tdat_mrg,1),1)]'))';
valtrans = valtrans(1:256, :);
diff     = vecnorm((valtrans(:, 1:3)-pgdat(1:256, :)), 2, 2);

% % % %--------------------------------------------------------------------------
% % % % Do an outlier detection
% % % outliers = find(diff>outthresh);
% % % if ~isempty(outliers)
% % %     outliers
% % %     diff(outliers)
% % %     error('stop')
% % % end
% % % tdat_mrg(outliers, :) = NaN;
% % % valtrans(outliers, :) = NaN;
good_els              = find(isnan(tdat_mrg(:, 1))==0);
n_el                  = length(good_els);

% Return the 
rms_with_filter       = sqrt(1/n_el*sum( diff(good_els).^2 ) );
mdn                   = median(abs(diff(good_els))); % Median error
max_err               = max(abs(diff(good_els)));


%--------------------------------------------------------------------------
if dbg_flg == 1
    %----------------------------------------------------------------------
    disp(['filtered rms= ', num2str(rms_with_filter*1000), ' mm'])
    disp(['number of electrodes = ', num2str(length(good_els))]);
    %save(['merging_data_', subj_name], 'mdn', 'rms_with_filter', 'max_err')
    % figure
    % hold on

    %----------------------------------------------------------------------
    % Plot the nominal eeg cap and the
    figure
    set(gcf,'position',[680         394        1109         584])
    subplot(1,2,1)
    hold on
    plot3(pgdat(:,1),pgdat(:,2), pgdat(:,3),'.k','markersize',16)
    plot3(valtrans(:,1),valtrans(:,2),valtrans(:,3),'.r','markersize',16)

    for n = 1:size(tdat_mrg,1)-1
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
%     azs = [0:45:(360-45)]
%     els = [15 30];
%     for n2 = 1:length(els)
%         for n1 = 1:length(azs)
%             set(gca,'view',[azs(n1) els(n2)]);
%
%         end
%     end


%     writerObj = VideoWriter('figs/mergefinal_pointcloud_mov.mp4','MPEG-4');
%     writerObj.FrameRate = 2;
%     open(writerObj);
%
%     azs   = (0:5:360)';
%     azels = [ ...
%         azs 22.5+0*(0:5:360)'; ...
%         azs 45+0*(0:5:360)'; ...
%         ];
%     for k = 1:size(azels,1)
%         for j = 1:2
%             subplot(1,2,j)
%             view(azels(k,:))
%         end
%         %-------------------------------------------------------------------
%         F = getframe(gcf);
%         writeVideo(writerObj,F);
%         pause(0.01)
%
%     end
%     close(writerObj);

tdat_filt=tdat_mrg;
end









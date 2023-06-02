function [tdat_mrg_fullset, rms, max_err] = find_missing_elec_labels_err_vf(nomdat,pgdat, tdat_mrg,colobj_mrg,dbg_flg, subj_name, file_ext)
%after all of the electrodes have been labeled, perform a best fit to 
% label the rest and test the error 
%very similar to calc_valerr 
tdat_mrg=tdat_mrg(1:256, :);
is_fnd = find(isnan(tdat_mrg(:,1))==0);

disp([num2str(length(is_fnd)),' Electrodes found'])

%-------------------------------------------------------------------------------
% PART 1: Fill in gaps in tdat_mrg using nominal cap
% Collect the target and corresponding nominal electrode points
A  = [tdat_mrg(is_fnd,:) ones(length(is_fnd),1)]';
B0  = [nomdat.eeg_mps(is_fnd,:)]';
%full set of nominal points 
C=[nomdat.eeg_mps, ones(size(nomdat.eeg_mps, 1),1)]';
%-------------------------------------------------------------------------------
%compute best fit transformation of tdat to nomdat
sclf = linspace(0.5,1.5,30);
[xscf,yscf,zscf]=meshgrid(sclf,sclf,sclf);
sclf_ps = [xscf(:) yscf(:) zscf(:)];
for k = 1:size(sclf_ps,1)
    Btmp = (diag(sclf_ps(k,:))*(B0));
        % B = [sclf(k)*B0 ones(length(is_fnd),1)]';
        B = [Btmp; ones(1, length(is_fnd))];
        %-----------------------------------------------------------------------
        % Get registered points from each frame
        [T, rmserr]   = transform(A, B, 1);
        Ts{k}         = T;
        rmserrs(k)    = rmserr;        
end

 % find best fitting matches
    [tmp,k]   = min(rmserrs);
%apply best fit transformation to fullset
    bestsclf = sclf_ps(k,:);

tdat_nominal_fullset=(Ts{k}*(C.*sclf_ps(k)))';
tdat_nominal_fullset=tdat_nominal_fullset(:, 1:3);
disp(['RMS error: ',num2str(rmserrs(k)*1000),' mm'])
%apply that same transformation to the whole nominal cap
tdat_nominal_fullset=tdat_nominal_fullset(:, 1:3);
%fill in tdat_mrg gaps with this nominal
%try walk labeling?
gaps = find(isnan(tdat_mrg(:,1))==1);
%project transformed nominal points onto the surface of the merged point cloud 
fill_in=tdat_nominal_fullset(gaps, :);
%calculate normal vectors of colobj closest to
% Set parameters
elrad = 5/1000; % radius
loc=colobj_mrg.Location;
%find closest point to each label
elthick=8/1000;
figure;
hold on
cent=mean(tdat_nominal_fullset);
% for k=1:size(fill_in, 1);
%     dist=vecnorm(fill_in(k, :)-loc, 2, 2);
%     [dist_sort, idx]=sort(dist);
%     closest_node=loc(idx(1), :);
%     %create unit vector that points from closest node to center
%     cent_vec=closest_node-cent./(vecnorm(closest_node-cent, 2, 2));
%     %multiply that by the electrode thickness to get new location
% new_loc(k, :)=closest_node-elthick*cent_vec;
%  plot3(new_loc(k, 1), new_loc(k, 2), new_loc(k, 3), 'b.', 'MarkerSize', 10);
%  hold on
%  plot3(closest_node(1), closest_node(2), closest_node(3), 'r.', 'MarkerSize', 10);
% end
tdat_mrg_fullset=tdat_mrg;
tdat_mrg_fullset(gaps, :)=fill_in;
%all of tdat is now filled in
%compute best fit transformation from tdat to pgdat again

if file_ext=='.sfp'
pgdat=pgdat.elpos/100;
else
pgdat=pgdat/100;
end
A_full  = [tdat_mrg_fullset(1:256, :) ones(256, 1)]';
B0_full  = [pgdat(1:256, :)]';
B_full = [B0_full; ones(1, 256)];
% for k = 1:size(sclf_ps,1)
%     Btmp = (diag(sclf_ps(k,:))*(B0_full));
%         % B = [sclf(k)*B0 ones(length(is_fnd),1)]';
%         B_full = [Btmp; ones(1, 256)];
%         %-----------------------------------------------------------------------
%         % Get registered points from each frame
[T, rmserr]   = transform(A_full, B_full, 1);
%         Ts_full{k}         = T;
%         rmserrs(k)    = rmserr;        
% end
 % find best fitting matches
    [tmp,k]   = min(rmserrs);
valtrans=(T*B_full)';
diff=vecnorm((valtrans(:, 1:3)-tdat_mrg(1:256, :)), 2, 2);

rms=mean(abs(diff(isnan(diff)==0)));
    mdn=median(abs(diff(isnan(diff)==0)));
    max_err=max(abs(diff));
    disp(['mean rms= ', num2str(rms*1000), ' mm']) 
    disp(['median rms= ', num2str(mdn*1000), ' mm']) 
    disp(['max rms= ', num2str(max_err*1000), ' mm']) 
    save(['fullset_merging_data_', subj_name], 'mdn', 'rms', 'max_err')
disp(['tdat full set / photogrammetry RMS Error: ',num2str(rmserrs(k)*1000),' mm'])
figure
hold on

%-------------------------------------------------------------------------------
if dbg_flg == 1
        
    %---------------------------------------------------------------------------
    % Plot the nominal eeg cap and the 
    figure
    set(gcf,'position',[680         394        1109         584])
    subplot(1,2,1)
    hold on
    plot3(tdat_mrg_fullset(:,1),tdat_mrg_fullset(:,2),tdat_mrg_fullset(:,3),'.k','markersize',16)
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
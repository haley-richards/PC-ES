%-------------------------------------------------------------------------------
% 
% Perform ICP across scans between unlabeled electrodes consecutively
% through the scans
% 
%-------------------------------------------------------------------------------
function allscans = cpd_merging(allscans,nomdat,dbg_flg)


%-------------------------------------------------------------------------------
% Merge everything into to the nominal frame
Nsc  = length(allscans);    
Nsmp = 1000;
if dbg_flg == 1
    figure
    set(gcf,'position',[ 355         290        1380         617])
end
for n = 1:Nsc
    %---------------------------------------------------------------------------
    % Perform ICP: We want to match the nominal points to the scanned points
    Tinit  = [eye(3) zeros(3,1); zeros(1,3) 1];
    Aptcld = pointCloud(allscans(n).elunlab);
    Bptcld = pointCloud(nomdat.eeg_mps);
    % A = double(allscans(n).elunlab');
    % B = double(nomdat.eeg_mps');
    % tform = pcregistercpd(moving,fixed);
    [tform,Bfinal] = pcregistercpd(Bptcld,Aptcld);
    % [Y,Tfinal,Dist_err,iter] = icp(A,B,Tinit);
    % disp(['dist_err from icp: ',num2str(1000*Dist_err),' mm'])
    % Tfinal
    %     whos
    %     norm(tform - Bfinal.Location)
    Bfinal = Bfinal.Location;
    %---------------------------------------------------------------------------
    % Use the labels from the ICP work, i.e. after the transformation match it
    % to the closest point. Since just finding the closest electrode may not 
    % result in a unique set. We force the labeled set to be unique, and
    % then randomly pick the assigning order with a bunch of samples. We
    % pick the best as the sample with the small overall RMS between sets
    % of points
    % fig
    % Bfinal  = (Tfinal*([nomdat.eeg_mps ones(size(nomdat.eeg_mps,1),1)]'))';
    Nelfnd  = size(allscans(n).elunlab,1);
    el_labs = zeros(Nelfnd,Nsmp);    
    rms_v   = zeros(Nsmp,1);
    for k1 = 1:Nsmp
        js   = randsample(Nelfnd,Nelfnd);
        is   = 1:size(Bfinal,1);
        for k2 = 1:Nelfnd            
            [ds,i]             = min( sum( (Bfinal(is,1:3)-repmat(allscans(n).elunlab(js(k2),:),length(is),1)).^2,2) );        
            el_labs(js(k2),k1) = is(i);
            is(i)              = [];
        end
        rms_v(k1) = mean(sqrt(sum( (Bfinal(el_labs(:,k1),1:3)-allscans(n).elunlab).^2,2)));
    end
    if dbg_flg == 1
        subplot(1,2,1);hold on
        plot(sort(rms_v))
        title('sorted(RMS) per sample')
        subplot(1,2,2);hold on
        plot(n,length(unique(rms_v))/length(rms_v)*100,'.','markersize',16)
        title('Percent Unique Labels')
    end
    
    %---------------------------------------------------------------------------
    % Pick out the labeling that minimizes the error
    [tmp,i] = min(rms_v);   
    el_labs = el_labs(:,i);      
    allscans(n).el_labs = el_labs;
end

% figure
% pcshow(colobj,'Markersize',MS);
% hold on
% elunlab = double(elunlab);
% plot3(elunlab(:,1),elunlab(:,2),elunlab(:,3),'.k','markersize',16)
% for n = 1:length(el_labs)
%     text(elunlab(n,1),elunlab(n,2),elunlab(n,3),num2str(el_labs(n)),'fontsize',12,'color','r')
% end


%-------------------------------------------------------------------------------
%
% Simultaneously merge all the scans into
%
%-------------------------------------------------------------------------------
function [allscans,nomdat] = bestfit_scan_to_nominal(allscans0,nomdat,dbg_flg)

%-------------------------------------------------------------------------------
% Set parameters
sclf = linspace(0.5,1.5,15);
[xscf,yscf,zscf]=meshgrid(sclf,sclf,sclf);
sclf_ps = [xscf(:) yscf(:) zscf(:)];
fac  = 1.05;
fac2 = 1.09;

%-------------------------------------------------------------------------------
Nsc     = length(allscans0);
Ts      = cell(size(sclf_ps,1),Nsc);
rmserrs = zeros(size(sclf_ps,1),Nsc);
n = 1;%
%---------------------------------------------------------------------------
% Get the indices of the located target points
tdat   = allscans0(n).tdat;
is_fnd = find( (isnan(tdat(:,1)) == 0) );
is_fnd = is_fnd( is_fnd < 257);

%---------------------------------------------------------------------------
% Collect the target and corresponding nominal electrode points
B  = [tdat(is_fnd,:) ones(length(is_fnd),1)]';
A0 = nomdat.eeg_mps(is_fnd,:);

%---------------------------------------------------------------------------
% For each scale factor find a best fit between the target points from
% the point cloud and the nominal electrodes points
for k = 1:size(sclf_ps,1)
    Atmp = (diag(sclf_ps(k,:))*(A0'))';
    A = [Atmp ones(length(is_fnd),1)]';
    %-----------------------------------------------------------------------
    % Get registered points from each frame
    % T is such that A = T * B
    [T, rmserr]   = transform(A, B, 1);
    Ts{k,n}       = T;
    rmserrs(k,n)  = rmserr;
    % rmserrs(k,n)  = max(max( A - (T*B)));
end


%-------------------------------------------------------------------------------
% Find the best fitting matches
rmserrsv  = sqrt(sum(rmserrs.^2,2));
[tmp,k]   = min(rmserrsv);
disp(['Best RMS of matching targets with nominal elecs = ',num2str(round(1000*rmserrs(k,:),1)),' mm'])
disp(['Best scale factor = ',num2str(sclf_ps(k,:))]);

%-------------------------------------------------------------------------------
% Scale the best fitting result
nomdat.eeg_mps = (diag(sclf_ps(k,:))*(nomdat.eeg_mps'))';

%-------------------------------------------------------------------------------
% Transform each scan into the nominal EEG frame
n = 1;
%---------------------------------------------------------------------------
tdat    = allscans0(n).tdat;
elunlab = allscans0(n).elunlab;
locs    = allscans0(n).colobj.Location;
T       = Ts{k,n};
%---------------------------------------------------------------------------
% Transform each thing: labeled electrodes, unlabeled electrodes, and the
% point cloud
B        = [tdat ones(size(tdat,1),1)]';
tdat     = (T*B)';
if ~isempty(elunlab)
    B        = [elunlab ones(size(elunlab,1),1)]';
    elunlab  = (T*B)';
end
B        = [locs ones(size(locs,1),1)]';
locs     = (T*B)';
colobj   = pointCloud(locs(:,1:3), 'Color', allscans0(n).colobj.Color);
%---------------------------------------------------------------------------
% Remove outliers from the unlabeled electrode list
if ~isempty(elunlab)
    ris = [];
    for i = 1:size(elunlab,1)
        ds = sqrt(sum( (nomdat.eeg_mps - repmat(elunlab(i,1:3),size(nomdat.eeg_mps,1),1)).^2,2));
        if min(ds) > 0.02
            ris = [ris; i];
        end
    end
    elunlab(ris,:) = [];
end

%---------------------------------------------------------------------------
% Record the transformed data
allscans(n).colobj  = colobj;
allscans(n).tdat    = tdat(:,1:3);
allscans(n).T       = T;
if ~isempty(elunlab)
    allscans(n).elunlab = elunlab(:,1:3);
else
    allscans(n).elunlab = [];
end
if isfield(allscans0(n),'tri') == 1
    allscans(n).tri = allscans0(n).tri;
end




%-------------------------------------------------------------------------------
if dbg_flg == 1
    figure
    plot(rmserrs)
    
    %---------------------------------------------------------------------------
    % Plot the nominal eeg cap and the
    figure
    set(gcf,'position',[680         394        1109         584])
    subplot(1,2,1)
    hold on
    plot3(nomdat.eeg_mps(:,1),nomdat.eeg_mps(:,2),nomdat.eeg_mps(:,3),'.k','markersize',12)
    for n = 1:Nsc
        elunlab = allscans(n).elunlab;
        if ~isempty(elunlab)
            plot3(elunlab(:,1),elunlab(:,2),elunlab(:,3),'.','markersize',16)
        end
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
    for n = 1:Nsc
        colobj = allscans(n).colobj;
        pcshow(colobj,'Markersize',MS);
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
    
    %     azs = [0:45:(360-45)]
    %     els = [15 30];
    %     for n2 = 1:length(els)
    %         for n1 = 1:length(azs)
    %             set(gca,'view',[azs(n1) els(n2)]);
    %
    %         end
    %     end
    return
    
    writerObj = VideoWriter('figs/merge_pointcloud_mov.mp4','MPEG-4');
    writerObj.FrameRate = 2;
    open(writerObj);
    
    azs   = (0:5:360)';
    azels = [ ...
        azs 22.5+0*(0:5:360)'; ...
        azs 45+0*(0:5:360)'; ...
        ];
    for k = 1:size(azels,1)
        for j = 1:2
            subplot(1,2,j)
            view(azels(k,:))
        end
        %-------------------------------------------------------------------
        F = getframe(gcf);
        writeVideo(writerObj,F);
        pause(0.01)
        
    end
    close(writerObj);
    
    
end







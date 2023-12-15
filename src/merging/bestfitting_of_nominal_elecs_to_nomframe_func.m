%-------------------------------------------------------------------------------
%
% Find a best fit between labeled electrodes from the scans and labeled
% electrodes from the nominal hdeeg locations. The nominal data is
% stretched and rotated and translated to fit the scan data.
%
%-------------------------------------------------------------------------------
function [tdat,nomdat,colobj,bestsclf,Top] = bestfitting_of_nominal_elecs_to_nomframe_func(tdat,nomdat,colobj,dbg_flg)

%-------------------------------------------------------------------------------
% Set parameters
sclf = linspace(0.5,1.5,15);
[xscf,yscf,zscf]=meshgrid(sclf,sclf,sclf);
sclf_ps = [xscf(:) yscf(:) zscf(:)];
fac  = 1.05;
fac2 = 1.09;

%-------------------------------------------------------------------------------
% Get the indices of the located target points
is_fnd = find( (isnan(tdat(:,1)) == 0) );
is_fnd = is_fnd( is_fnd < 257);


%---------------------------------------------------------------------------
% Collect the target and corresponding nominal electrode points
A  = [tdat(is_fnd,:) ones(length(is_fnd),1)]';
B0 = nomdat.eeg_mps(is_fnd,:);

%---------------------------------------------------------------------------
% For each scale factor find a best fit between the target points from
% the point cloud and the nominal electrodes points
rmserrs = zeros(size(sclf_ps,1),1);
for k = 1:size(sclf_ps,1)
    Btmp = (diag(sclf_ps(k,:))*(B0'))';
    % B = [sclf(k)*B0 ones(length(is_fnd),1)]';
    B = [Btmp ones(length(is_fnd),1)]';
    %-----------------------------------------------------------------------
    % Get registered points from each frame
    [T, rmserr]   = transform_loc(B, A, 1);
    Ts{k}         = T;
    rmserrs(k)    = rmserr;
end
%     figure
%     plot(rmserrs)
%
%-------------------------------------------------------------------------------
% Transform to the best-fit space
%   1. stretch the nominal to the best fitting one
%   2. transform the target data and point cloud to the nominal frame
[tmp,k]  = min(rmserrs);
Top      = Ts{k};
bestsclf = sclf_ps(k,:);
disp(['Best RMS of matching targets with nominal elecs = ',num2str(round(1000*rmserrs(k),1)),' mm', ...
    ', best scale factor = ',num2str(sclf_ps(k,:))]);          % ', best scale factor = ',num2str(sclf(k))]);
%   1. stretch the nominal to the best fitting one
nomdat.eeg_mps = (diag(sclf_ps(k,:))*(nomdat.eeg_mps'))';
%   2. transform the target data and point cloud to the nominal frame
locs   = colobj.Location;
tdat   = (Top*([tdat ones(size(tdat,1),1)]'))';
locs   = (Top*([locs ones(size(locs,1),1)]'))';
tdat   = tdat(:,1:3);
colobj = pointCloud(locs(:,1:3), 'Color', colobj.Color);



%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
if dbg_flg == 1
    figure
    hold on
    pcshow(colobj)
    plot3(nomdat.eeg_mps(:,1),nomdat.eeg_mps(:,2),nomdat.eeg_mps(:,3),'.c','markersize',16)
    plot3(tdat(:,1),tdat(:,2),tdat(:,3),'.r','markersize',16)
    for n = 1:length(is_fnd)
        text(fac2*tdat(is_fnd(n),1),fac2*tdat(is_fnd(n),2),fac2*tdat(is_fnd(n),3),num2str(is_fnd(n)),'fontsize',12,'color','c')
    end
end




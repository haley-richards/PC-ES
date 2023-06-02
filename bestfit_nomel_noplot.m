%--------------------------------------------------------------------------
%
%
%
%--------------------------------------------------------------------------
function [Bout,elns,cnxs,bestsclf] = bestfit_nomel_noplot(tdat,nomdat,colobj,sclf_ps)

%--------------------------------------------------------------------------
% Set parameters
if isempty(sclf_ps)
    sclf = linspace(0.5,1.5,15);
    [xscf,yscf,zscf]=meshgrid(sclf,sclf,sclf);
    sclf_ps = [xscf(:) yscf(:) zscf(:)];
end

%-------------------------------------------------------------------------------
% Get the indices of the located target points
is_fnd = find( (isnan(tdat(:,1)) == 0) );
is_fnd = is_fnd( is_fnd < 257);


% %---------------------------------------------------------------------------
% % Collect the target and corresponding nominal electrode points
% A  = [tdat(is_fnd,:) ones(length(is_fnd),1)]';
% B0 = nomdat.eeg_mps(is_fnd,:);
% 
% %---------------------------------------------------------------------------
% % For each scale factor find a best fit between the target points from
% % the point cloud and the nominal electrodes points
% Btmp = (diag(sclf_ps)*(B0'))';
% % B = [sclf(k)*B0 ones(length(is_fnd),1)]';
% B = [Btmp ones(length(is_fnd),1)]';
% %-----------------------------------------------------------------------
% % Get registered points from each frame
% T    = transform(A, B, 1);
% Btmp = (diag(sclf_ps)*(nomdat.eeg_mps'))';
% Bout = double((T*([Btmp ones(size(nomdat.eeg_mps,1),1)]'))');
% elns = 1:257;

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
    [T, rmserr]   = transform(A, B, 1);
    Ts{k}         = T;
    rmserrs(k)    = rmserr;
end

%---------------------------------------------------------------------------
% Plot the best fitting matches
[tmp,k]   = min(rmserrs);
% set(handles.message_texts,'String',['Best RMS of matching targets with nominal elecs = ',num2str(round(1000*rmserrs(k),1)),' mm', ...
%     ', best scale factor = ',num2str(sclf_ps(k,:))]);          % ', best scale factor = ',num2str(sclf(k))]);
if size(sclf_ps,1) > 3
 disp(['Best RMS of matching targets with nominal elecs = ',num2str(round(1000*rmserrs(k),1)),' mm', ...
     ', best scale factor = ',num2str(sclf_ps(k,:))]);
end
Btmp     = (diag(sclf_ps(k,:))*(nomdat.eeg_mps'))';
Bout     = double((Ts{k}*([Btmp ones(size(nomdat.eeg_mps,1),1)]'))');
% Btmp     = (diag(sclf_ps(k,:))*(handles.eeg_mps'))';
% Bout     = double((Ts{k}*([Btmp ones(size(handles.eeg_mps,1),1)]'))');
elns     = 1:258;
bestsclf = sclf_ps(k,:);

%---------------------------------------------------------------------------
% Find the electrodes that are within the bounds of the point cloud
dxyz = 45;
is_keep = find( ...
    (min(colobj.Location(:,1))-dxyz/1000 < Bout(:,1)) & (Bout(:,1) < max(colobj.Location(:,1))+dxyz/1000 ) & ...
    (min(colobj.Location(:,2))-dxyz/1000 < Bout(:,2)) & (Bout(:,2) < max(colobj.Location(:,2))+dxyz/1000 ) & ...
    (min(colobj.Location(:,3))-dxyz/1000 < Bout(:,3)) & (Bout(:,3) < max(colobj.Location(:,3))+dxyz/1000 ) );
Bout = Bout(is_keep,:);
elns = elns(is_keep);
ris  = [];
for n = 1:length(elns)
    ds = sqrt( sum( (double(colobj.Location) - repmat(Bout(n,1:3),size(colobj.Location,1),1)).^2,2));
    if min(ds)*1000 > 30
        ris = [ris; n];
    end
end
Bout(ris,:) = [];
elns(ris)   = [];

%---------------------------------------------------------------------------
% Keep only the connection with the remaining electrodes
cnxs = nomdat.cnxs;
kis  = []; 
for n = 1:size(cnxs,1)
    if ~isempty(intersect(cnxs(n,1),elns)) && ~isempty(intersect(cnxs(n,2),elns)) 
        kis  = [kis; n];
    end
end
cnxs = cnxs(kis,:);
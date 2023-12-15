%-------------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------------
function [Bout,T,best_scf] = bestfit_tdat_to_nomel_noplot_fullset_loc(tdat,nomdat,dbg_flg)

%--------------------------------------------------------------------------
% Set parameters
sclf = linspace(0.5,1.5,15);
[xscf,yscf,zscf]=meshgrid(sclf,sclf,sclf);
sclf_ps = [xscf(:) yscf(:) zscf(:)];

%--------------------------------------------------------------------------
% Get the indices of the located target points
is_fnd = find( (isnan(tdat(:,1)) == 0) );
is_fnd = is_fnd( is_fnd < 257);


%--------------------------------------------------------------------------
% Collect the target and corresponding nominal electrode points
B  = [tdat(is_fnd,:) ones(length(is_fnd),1)]';
A0 = nomdat.eeg_mps(is_fnd,:);

%--------------------------------------------------------------------------
% For each scale factor find a best fit between the target points from
% the point cloud and the nominal electrodes points
rmserrs = zeros(size(sclf_ps,1),1);
for k = 1:size(sclf_ps,1)
    Atmp = (diag(sclf_ps(k,:))*(A0'))';
    % B = [sclf(k)*B0 ones(length(is_fnd),1)]';
    A = [Atmp ones(length(is_fnd),1)]';
    %----------------------------------------------------------------------
    % Get registered points from each frame
    [T, rmserr]   = transform_loc(A, B, 1);
    Ts{k}         = T;
    rmserrs(k)    = rmserr;
end

%--------------------------------------------------------------------------
% Plot the best fitting matches
[tmp,k]   = min(rmserrs);
% set(handles.message_texts,'String',['Best RMS of matching targets with nominal elecs = ',num2str(round(1000*rmserrs(k),1)),' mm', ...
%     ', best scale factor = ',num2str(sclf_ps(k,:))]);          % ', best scale factor = ',num2str(sclf(k))]);
if dbg_flg == 1
    if size(sclf_ps,1) > 3
        disp(['Best RMS of matching targets with nominal elecs = ',num2str(round(1000*rmserrs(k),1)),' mm', ...
            ', best scale factor = ',num2str(sclf_ps(k,:))]);
    end
end
best_scf = sclf_ps(k,:);
% Btmp     = (diag(sclf_ps(k,:))*(nomdat.eeg_mps'))';
Bout     = double((Ts{k}*([tdat ones(size(tdat,1),1)]'))');
Bout     = Bout(:,1:3);
T        = Ts{k};
%--------------------------------------------------------------------------
bestsclf = sclf_ps(k,:);


%-------------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------------
function [Bout,elns,bestsclf] = try_bestfitting_of_nominal_elecs(handles)

%-------------------------------------------------------------------------------
% Set parameters
if isfield(handles,'bestsclf') == 1
    sclf_ps = handles.bestsclf;
else
    sclf = linspace(0.5,1.5,15);
    [xscf,yscf,zscf]=meshgrid(sclf,sclf,sclf);
    sclf_ps = [xscf(:) yscf(:) zscf(:)];
end
fac  = 1.05;
fac2 = 1.09;
%-------------------------------------------------------------------------------
% Perform translations of the points
tic
[colobj,tdat,eclks] = update_PLY_pts(handles);
toc

%-------------------------------------------------------------------------------
% Get the indices of the located target points
is_fnd = find( (isnan(tdat(:,1)) == 0) );
is_fnd = is_fnd( is_fnd < 258);

if length(is_fnd) < 3
    set(handles.message_texts,'String','Not enough target points found. Find some more!')    
else
    %---------------------------------------------------------------------------
    % Collect the target and corresponding nominal electrode points
    A  = [tdat(is_fnd,:) ones(length(is_fnd),1)]';    
    B0 = handles.eeg_mps(is_fnd,:);  
    
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
    %     figure
    %     plot(rmserrs)
    %
    %---------------------------------------------------------------------------
    % Plot the best fitting matches
    [tmp,k]   = min(rmserrs);
    set(handles.message_texts,'String',['Best RMS of matching targets with nominal elecs = ',num2str(round(1000*rmserrs(k),1)),' mm', ...
        ', best scale factor = ',num2str(sclf_ps(k,:))]);          % ', best scale factor = ',num2str(sclf(k))]);        
    Btmp = (diag(sclf_ps(k,:))*(handles.eeg_mps'))';    
    % Bout = double((Ts{k}*([ sclf(k)*handles.eeg_mps ones(size(handles.eeg_mps,1),1)]'))');
    Bout = double((Ts{k}*([Btmp ones(size(handles.eeg_mps,1),1)]'))');
    elns = 1:257;
    bestsclf = sclf_ps(k,:);
    
    %---------------------------------------------------------------------------
    % Find the electrodes that are within the bounds of the point cloud
    is_keep = find( ...
        (min(colobj.Location(:,1)) < Bout(:,1)) & (Bout(:,1) < max(colobj.Location(:,1)) ) & ... 
        (min(colobj.Location(:,2)) < Bout(:,2)) & (Bout(:,2) < max(colobj.Location(:,2)) ) & ...
        (min(colobj.Location(:,3)) < Bout(:,3)) & (Bout(:,3) < max(colobj.Location(:,3)) ) );
    Bout = Bout(is_keep,:);
    elns = elns(is_keep);    
    ris  = [];
    for n = 1:length(elns)
        ds = sqrt( sum( (double(colobj.Location) - repmat(Bout(n,1:3),size(colobj.Location,1),1)).^2,2));
        if min(ds)*1000 > 10
            ris = [ris; n];
        end
    end
    Bout(ris,:) = [];
    elns(ris)   = [];
    
    %---------------------------------------------------------------------------
    axes(handles.mainax);
    hold on
    plot3(fac*Bout(:,1),fac*Bout(:,2),fac*Bout(:,3),'.c','markersize',16)
    for n = 1:length(elns)
        text(fac2*Bout(n,1),fac2*Bout(n,2),fac2*Bout(n,3),num2str(elns(n)),'fontsize',12,'color','c')
    end
end


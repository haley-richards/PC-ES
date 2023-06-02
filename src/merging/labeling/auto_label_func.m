%--------------------------------------------------------------------------
%
% Auto labeling algorithm.
%
%--------------------------------------------------------------------------
function [tdat,n_unlab,n_ellab,elunlab] = auto_label_func(elunlab,tdat,colobj,nomdat,dethr,dbg_flg)

%--------------------------------------------------------------------------
% dethr   = 8/1000; % Maximum distance to consider for automatic labeling of
% electrodes (m)


%--------------------------------------------------------------------------
% 1. Perform a best-fit between the nominal and registration data.
nomps0 = bestfit_nomel_noplot_fullset(tdat,nomdat);
if dbg_flg == 1
    update_resplot(elunlab,tdat,colobj,nomps0)
    title([num2str(size(elunlab,1)),' unlabeled, ',num2str(length(find( (isnan(tdat(:,1)) == 0) ))),' labeled electrodes'])
end


%--------------------------------------------------------------------------
check = 1;
while check == 1
    check = 0;
    %--------------------------------------------------------------------------
    % 2. Label the electrodes associated with the closest fitting electrodes
    is_fnd    = find( (isnan(tdat(:,1)) == 0) );
    is_fnd    = is_fnd( is_fnd < 257);
    ellab     = tdat(is_fnd,:);
    nomps     = nomps0(is_fnd,:);
    dists     = sqrt(sum( ellab - nomps,2).^2);
    [tmp,sid] = sort(dists,'ascend');
    elchcks   = is_fnd(sid);
    % Loop through all the electrodes to be checked
    for n = 1:length(elchcks)
        if dbg_flg == 1        
            disp(['Working on E',num2str(elchcks(n)),', ',num2str(length(is_fnd)),' El. labeled'])
        end
        %----------------------------------------------------------------------
        % If the current electrode is less than the threshold for labeling,
        % then label adjacent electrodes if they are close enough to nominal
        % electrodes
        ds      = sqrt(sum( (nomps0 - repmat(tdat(elchcks(n),:),size(nomps0,1),1)).^2,2));
        [mind,i] = min(ds);
        if mind < dethr
            % Get the unlabeled neighbors of current electrode
            is     = find( (nomdat.cnxs(:,1) == elchcks(n)) | (nomdat.cnxs(:,2) == elchcks(n)) );
            adj_es = setdiff( unique( [nomdat.cnxs(is,1);nomdat.cnxs(is,2)]), elchcks );
            for k = 1:length(adj_es)
                % disp(['   Adjacent E',num2str(adj_es(k)),' of ',num2str(length(adj_es))])
                
                % Label the electrode with the closest nominal electrode
                ds      = sqrt(sum( (elunlab - repmat(nomps0(adj_es(k),:),size(elunlab,1),1)).^2,2));
                [mind,i] = min(ds);
                if mind < dethr
                    tdat(adj_es(k),:) =  elunlab(i,:);
                    elunlab(i,:)      = [];
                    check = 1;
                end
            end
        end
    end
    %------------------------------------------------------------------
    % 1. Perform a best-fit between the nominal and registration data.
    nomps0 = bestfit_nomel_noplot_fullset(tdat,nomdat);
   
    
end

if dbg_flg == 1
    update_resplot(elunlab,tdat,colobj,nomps0)
    title([num2str(size(elunlab,1)),' unlabeled, ',num2str(length(find( (isnan(tdat(:,1)) == 0) ))),' labeled electrodes'])
end
n_unlab = size(elunlab,1);
n_ellab = length(find( (isnan(tdat(:,1)) == 0)));

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Subfunctions
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Update the results
%
%--------------------------------------------------------------------------
function update_resplot(elunlab,tdat,colobj,nomps)
MS = 200; % Marksize of point cloud
figure
pcshow(colobj,'Markersize',MS);
hold on
%--------------------------------------------------------------------------
% Nominal electrodes
fac  = 1.05;
fac2 = 1.09;
plot3(fac*nomps(:,1),fac*nomps(:,2),fac*nomps(:,3),'.c','markersize',16)
for n = 1:size(nomps,1)
    text(fac2*nomps(n,1),fac2*nomps(n,2),fac2*nomps(n,3),num2str(n),'fontsize',12,'color','c')
end
%--------------------------------------------------------------------------
% Unlabeled electrodes
plot3(elunlab(:,1),elunlab(:,2),elunlab(:,3),'.m','markersize',32)
%--------------------------------------------------------------------------
% Labeled electrodes
plot3(tdat(:,1),tdat(:,2),tdat(:,3),'.y','markersize',36)
for n = 1:size(tdat,1)
    if isnan(tdat(n,1)) == 0
        tmp = double(tdat(n,:));
        text(1.02*tmp(1),1.02*tmp(2),1.02*tmp(3),num2str(n),'fontsize',16,'FontName','times','color','red')
    end
end
%--------------------------------------------------------------------------
% Formatting
axis equal
grid on
box on
FS = 12;
xlabel('X','fontsize',FS,'FontName','times')
ylabel('Y','fontsize',FS,'FontName','times')
zlabel('Z','fontsize',FS,'FontName','times')
set(gca,'FontSize',FS,'FontName','times')

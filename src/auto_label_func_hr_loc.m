%--------------------------------------------------------------------------
%
% Auto labeling algorithm.
%   Inputs:
%       elunlab - unlabeled electrodes
%          tdat - labeled electrodes
%        colobj - point clound
%        nomdat - nominal data
%         dethr - distance threshold (meters, same units as point cloud)
%   Outputs:
%          tdat - labeled electrodes (there should be a lot more)
%       n_unlab - number of unlabeled electrodes
%       n_ellab - number of labeled electrodes
%       elunlab - unlabeled electrodes (should be smaller)
% 
% Algorithm steps:
%   1. Perform a best-fit between the nominal and registration data.
%   2. Start the while loop, continuing to look for electrodes to label
%      as long as in each main loop there was at least one electrode newly
%      labeled.
%   2.b. Get a sorted list of indices that best match between     
%        the iPhone scan and nominal electrodes
%   2.c. Loop through all the electrodes to be checked
%   2.d. If the current electrode is less than the threshold for     
%        labeling, then label adjacent electrodes if they are close      
%        enough to nominal electrodes
%   2.d.i. Get the unlabeled neighbors of current electrode - using the
%          neighbors defined by the nominal mesh edges
%   2.d.ii.If there is an unlabeled electrode within the                
%          distance threshold of the kth adjacent electrode,                 
%          then label it accordingly and remove it from the                 
%          unlabeled electrode list. 
%   2.e. Update the best-fit between the nominal and registration data.
% 
%--------------------------------------------------------------------------
function [tdat,n_unlab,n_ellab,elunlab] = auto_label_func_hr_loc(elunlab,tdat,colobj,nomdat,dethr,dbg_flg,tri)

%--------------------------------------------------------------------------
% 1. Perform a best-fit between the nominal and registration data.
nomps0 = bestfit_nomel_noplot_fullset_loc(tdat,nomdat,dbg_flg);
if dbg_flg == 1
    update_resplot(elunlab,tdat,colobj,nomps0)
    title([num2str(size(elunlab,1)),' unlabeled, ',num2str(length(find( (isnan(tdat(:,1)) == 0) ))),' labeled electrodes'])
end


%--------------------------------------------------------------------------
%   2. Start the while loop, continuing to look for electrodes to label
%      until we hit a stopping criteria
check = 1;
iter  = 0;
while check == 1
    check = 0;
    %----------------------------------------------------------------------
    % 2.b. Get a sorted list of indices that best match between 
    %      the iPhone scan and nominal electrodes
    is_fnd    = find( (isnan(tdat(:,1)) == 0) ); % Get indices of found electrodes
    is_fnd    = is_fnd( is_fnd < 257);
    ellab     = tdat(is_fnd,:);                 % Labeled positions
    nomps     = nomps0(is_fnd,:);
    dists     = sqrt(sum( ellab - nomps,2).^2); % Distance between nominal and labeled
    [tmp,sid] = sort(dists,'ascend');
    elchcks   = is_fnd(sid);                    % ordered set of electrodes to check (closest to farthest)

    %----------------------------------------------------------------------
    % 2.c. Loop through all the electrodes to be checked
    for n = 1:length(elchcks)
        if dbg_flg == 1        
            disp(['Working on E',num2str(elchcks(n)),', ',num2str(length(is_fnd)),' El. labeled'])
        end
        %------------------------------------------------------------------
        % 2.d. If the current electrode is less than the threshold for 
        %      labeling, then label adjacent electrodes if they are close 
        %      enough to nominal electrodes
        ds      = sqrt(sum( (nomps0 - repmat(tdat(elchcks(n),:),size(nomps0,1),1)).^2,2));
        [mind,icls] = min(ds);
        if mind < dethr
            %--------------------------------------------------------------
            % 2.d.i. Get the unlabeled neighbors of current electrode
            is     = find( (nomdat.cnxs(:,1) == elchcks(n)) | (nomdat.cnxs(:,2) == elchcks(n)) );
            adj_es = setdiff( unique( [nomdat.cnxs(is,1);nomdat.cnxs(is,2)]), elchcks );
            for k = 1:length(adj_es)
                % disp(['   Adjacent E',num2str(adj_es(k)),' of ',num2str(length(adj_es))])
                % 2.d.ii. If there is an unlabeled electrode within the
                %         distance threshold of the kth adjacent electrode, 
                %         then label it accordingly and remove it from the 
                %         unlabeled electrode list. 
                ds      = sqrt(sum( (elunlab - repmat(nomps0(adj_es(k),:),size(elunlab,1),1)).^2,2));
                [mind,i] = min(ds);
                if mind < dethr
                    tdat(adj_es(k),:) =  elunlab(i,:);
                    elunlab(i,:)      = [];
                    check = 1;  % Continue the loop as long as at least one electrode was labeled

                    %------------------------------------------------------
                    % Plot what's going on. 
                    if dbg_flg == 1
                        nomdat.sens_ps = nomps0;                        
                        [xm,ym,zm] = sphere(20);                        
                        %--------------------------------------------------
                        % Illustrating the fit
                        % * Plot scan
                        % * Local nominal mesh
                        % * Electrodes labeled
                        % * Matching nominal electrodes
                        figure;hold on
                        set(gcf,'position',[1130         395        1174         826])
                        %--------------------------------------------------
                        % Plot the point cloud as a surface that is
                        % somewhat transparent
                        locs     = double(colobj.Location);
                        col_nds  = double(colobj.Color)/255;
                        col_face = zeros(1,size(tri,1),3);
                        col_face(1,:,1) = 1/3*( col_nds(tri(:,1),1)+col_nds(tri(:,2),1)+col_nds(tri(:,3),1) );
                        col_face(1,:,2) = 1/3*( col_nds(tri(:,1),2)+col_nds(tri(:,2),2)+col_nds(tri(:,3),2) );
                        col_face(1,:,3) = 1/3*( col_nds(tri(:,1),3)+col_nds(tri(:,2),3)+col_nds(tri(:,3),3) );
                        patch( ...
                            [locs(tri(:,1),1) locs(tri(:,2),1) locs(tri(:,3),1)]', ...
                            [locs(tri(:,1),2) locs(tri(:,2),2) locs(tri(:,3),2)]', ...
                            [locs(tri(:,1),3) locs(tri(:,2),3) locs(tri(:,3),3)]', ...
                            col_face,'linestyle','none','facealpha',0.6)                        
                        for i = 1:length(elchcks)                        
                            plot_docbrown_dstthresh(nomdat,'b',tdat(elchcks(i),:),0.09)
                        end
                        %--------------------------------------------------
                        plot3(tdat(elchcks,1),tdat(elchcks,2),tdat(elchcks,3),'.g','markersize',12)
                        plot3(nomps0(elchcks,1),nomps0(elchcks,2),nomps0(elchcks,3),'om','markersize',8,'linewidth',2)
                        %--------------------------------------------------
                        % Plot the threshold around each
                        for i = 1:length(elchcks)
                            surf(xm*dethr + tdat(elchcks(i),1),ym*dethr + tdat(elchcks(i),2),zm*dethr + tdat(elchcks(i),3), ...
                                'linestyle','none','facecolor','yellow','facealpha',0.5)
                        end

                        %--------------------------------------------------
                        camlight(0, 70)
                        view([70 42])
                        saveas(gcf,'figs/autolab_ill_registration','png')


                        %--------------------------------------------------
                        % Illustrating the algorithm
                        % * Plot scan
                        % * Local nominal mesh
                        % * Plot the current labeled electrode (from scan)
                        % * Plot the matching nominal electrode
                        % * Plot the adjacent edges
                        figure;hold on
                        set(gcf,'position',[1130         395        1174         826])
                        %--------------------------------------------------
                        % Plot the point cloud as a surface that is
                        % somewhat transparent
                        locs     = double(colobj.Location);
                        col_nds  = double(colobj.Color)/255;
                        col_face = zeros(1,size(tri,1),3);
                        col_face(1,:,1) = 1/3*( col_nds(tri(:,1),1)+col_nds(tri(:,2),1)+col_nds(tri(:,3),1) );
                        col_face(1,:,2) = 1/3*( col_nds(tri(:,1),2)+col_nds(tri(:,2),2)+col_nds(tri(:,3),2) );
                        col_face(1,:,3) = 1/3*( col_nds(tri(:,1),3)+col_nds(tri(:,2),3)+col_nds(tri(:,3),3) );
                        patch( ...
                            [locs(tri(:,1),1) locs(tri(:,2),1) locs(tri(:,3),1)]', ...
                            [locs(tri(:,1),2) locs(tri(:,2),2) locs(tri(:,3),2)]', ...
                            [locs(tri(:,1),3) locs(tri(:,2),3) locs(tri(:,3),3)]', ...
                            col_face,'linestyle','none','facealpha',0.6)                        
                        %--------------------------------------------------       
                        % Plot the edge connections
                        plot_docbrown_dstthresh(nomdat,'b',tdat(elchcks(n),:),0.09)
                        %--------------------------------------------------       
                        % Plot the current node of interest
                        plot3(tdat(elchcks(n),1),tdat(elchcks(n),2),tdat(elchcks(n),3),'.g','markersize',20)
                        plot3(tdat(elchcks,1),tdat(elchcks,2),tdat(elchcks,3),'.g','markersize',12)
                        %--------------------------------------------------
                        % Plot the closest associated nominal node and
                        % adjacent ones
                        plot3(nomps0(icls,1),nomps0(icls,2),nomps0(icls,3),'om','markersize',8,'linewidth',2)
                        plot3(nomps0(adj_es,1),nomps0(adj_es,2),nomps0(adj_es,3),'oc','markersize',8,'linewidth',2)
                        %--------------------------------------------------
                        % Plot the connected edges
                        plot3( ...
                            [nomps0(nomdat.cnxs(is,1),1) nomps0(nomdat.cnxs(is,2),1)]', ...
                            [nomps0(nomdat.cnxs(is,1),2) nomps0(nomdat.cnxs(is,2),2)]', ...
                            [nomps0(nomdat.cnxs(is,1),3) nomps0(nomdat.cnxs(is,2),3)]', ...
                            '-m','linewidth',2)
                        %--------------------------------------------------
                        % Plot the threshold
                        surf(xm*dethr + nomps0(adj_es(k),1),ym*dethr + nomps0(adj_es(k),2),zm*dethr + nomps0(adj_es(k),3), ...
                            'linestyle','none','facecolor','yellow','facealpha',0.5)

                        axis([tdat(elchcks(n),1)+60/1000*[-1 1] ...
                            tdat(elchcks(n),2)+60/1000*[-1 1] ...
                            tdat(elchcks(n),3)+60/1000*[-1 1] ])
                        camlight(0, 70)
                        view([70 42])
                        saveas(gcf,'figs/autolab_ill_adjacent_zoomin','png')

                        clear
                    end
                end
            end
        end
    end
    %----------------------------------------------------------------------
    % 2.e. Update the best-fit between the nominal and registration data.
    nomps0 = bestfit_nomel_noplot_fullset_loc(tdat,nomdat,dbg_flg);
   
    iter  = iter +1 ;
    
end

if dbg_flg == 1
    disp(['Number of iterations: ',num2str(iter)])
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

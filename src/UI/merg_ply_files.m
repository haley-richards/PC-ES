%-------------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------------
function [colobj_mrgs, tdat_mrg,nconxs,rmserrs,id_mrgs] = merg_ply_files(colobjs,tdats,eclkss)

MS = 300;
%-------------------------------------------------------------------------------
if nargin < 3
    for n1 = 1:length(colobjs)
        eclkss{n1} = [];
    end
elseif nargin < 2
    error('target information is needed')
end

%-------------------------------------------------------------------------------
% Sequentially combine the point clouds
nt       = size(tdats{1},1);
tdat_mrg = NaN*ones(nt,3);
ndat_mrg =    zeros(nt,1);
nconxs   = zeros(length(colobjs)-1,1);   
id_scans = 1:length(colobjs);
for n = 1:(length(colobjs)-1)   
    %---------------------------------------------------------------------------
    if n == 1
        %-----------------------------------------------------------------------
        % In the initial step combine the scans with the largest number of
        % intersecting targets/stickers
        %-----------------------------------------------------------------------
        % construct the connection matrix
        cnxs = zeros(length(colobjs));
        for n1 = 1:(length(colobjs)-1)
            for n2 = (n1+1):length(colobjs)
                is1 = find( isnan(tdats{n1}(:,1)) == 0);
                is2 = find( isnan(tdats{n2}(:,1)) == 0);
                cnxs(n1,n2) = length(intersect(is1,is2));
            end
        end
        %-----------------------------------------------------------------------
        % Get the location of the maximum combination
        is_max = find( max(cnxs(:)) == cnxs(:));
        [I,J]  = ind2sub(size(cnxs),is_max(1));
        nconxs(n) = cnxs(I,J);
        %-----------------------------------------------------------------------
        A   = tdats{I};
        B   = tdats{J};
    else
        %-----------------------------------------------------------------------
        % In subsequent iterations, find the largest number of intersecting
        % targets/stickers with the merged scan
        %-----------------------------------------------------------------------
        % construct the connection matrix
        cnxs = zeros(length(colobjs),1);
        is2  = find( isnan(tdat_mrg(:,1)) == 0);
        for n1 = 1:length(colobjs)
                is1      = find( isnan(tdats{n1}(:,1)) == 0);                
                cnxs(n1) = length(intersect(is1,is2));
        end        
        [tmp,I] = max(cnxs);
        nconxs(n) = tmp;
        %-----------------------------------------------------------------------
        A   = tdat_mrg;
        B   = tdats{I};
    end
    
    %---------------------------------------------------------------------------
    %---------------------------------------------------------------------------
    % Get registered points from each frame
    ikp           = find( (isnan(A(:,1)) == 0) & (isnan(B(:,1)) == 0) );
    A             = [A(ikp,:) ones(length(ikp),1)]';
    B             = [B(ikp,:) ones(length(ikp),1)]';
    [T, rmserr]   = transform_loc(A, B, 1);
    Bout          = T*B;
    rmserrs(n)    = rmserr;
    num_common(n) = length(ikp);
    %     figure;hold on
    %     pcshow(colobjs{I},'Markersize',MS);
    %     plot3(A(1,:),A(2,:),A(3,:),'.r','markersize',16)
    %     plot3(Bout(1,:),Bout(2,:),Bout(3,:),'.b','markersize',16)
    %     legend('A','B')
    %---------------------------------------------------------------------------
    
    %---------------------------------------------------------------------------
    % Merge target points, electrodes and point clounds
    %---------------------------------------------------------------------------        
    if n == 1
        tdat_out = (T*([tdats{J} ones(size(tdats{J},1),1)]'))';    
        for k = 1:nt
            if (isnan(tdats{I}(k,1)) == 0) && (isnan(tdat_out(k,1)) == 0) 
                tdat_mrg(k,:) = 1/2*(tdats{I}(k,:) + tdat_out(k,1:3));
                ndat_mrg(k)   = 2;
            elseif (isnan(tdats{I}(k,1)) == 0) 
                tdat_mrg(k,:) = tdats{I}(k,:);
                ndat_mrg(k)   = 1;
            elseif (isnan(tdat_out(k,1)) == 0)                  
                tdat_mrg(k,:) = tdat_out(k,1:3);
                ndat_mrg(k)   = 1;
            end
        end
    else
        tdat_out = (T*([tdats{I} ones(size(tdats{I},1),1)]'))';    
        for k = 1:nt
            if (ndat_mrg(k) > 0) && (isnan(tdat_out(k,1)) == 0) 
                tdat_mrg(k,:) = (ndat_mrg(k)*tdat_mrg(k,:) + tdat_out(k,1:3))/(ndat_mrg(k) + 1);
                ndat_mrg(k)   = ndat_mrg(k) + 1;
            elseif isnan(tdat_out(k,1)) == 0
                tdat_mrg(k,:) = tdat_out(k,1:3);
                ndat_mrg(k)   = 1;
            end
        end
    end
    % Merge electrodes
    if nargin == 3
        if n == 1
            eclks_out  = (T*([eclkss{J} ones(size(eclkss{J},1),1)]'))';
            eclks_mrg  = [eclkss{I}; eclks_out(:,1:3)];
        else
            eclks_out  = (T*([eclkss{I} ones(size(eclkss{I},1),1)]'))';
            eclks_mrg  = [eclks_mrg; eclks_out(:,1:3)];
        end
    end
    % Merge point clouds?!
    if n == 1
        locs_out   = (T*([colobjs{J}.Location ones(size(colobjs{J}.Location,1),1)]'))';
        locs_mrg   = [colobjs{I}.Location; locs_out(:,1:3)];
        col_mrg    = [colobjs{I}.Color; colobjs{J}.Color];
    else
        locs_out   = (T*([colobjs{I}.Location ones(size(colobjs{I}.Location,1),1)]'))';
        locs_mrg   = [locs_mrg; locs_out(:,1:3)];
        col_mrg    = [col_mrg; colobjs{I}.Color];
    end
    colobj_mrgs{n} = pointCloud(locs_mrg, 'Color', col_mrg);
    
    %---------------------------------------------------------------------------
    % Remove the individual data that was merged
    
    if n == 1
        id_mrgs{n}     = [I J];
        colobjs([I J]) = [];
        eclkss([I J])  = [];
        tdats([I J])   = [];
        id_scans([I J])= [];
    else
        id_mrgs{n} = id_scans(I);
        colobjs(I) = [];
        eclkss(I)  = [];
        tdats(I)   = [];
        id_scans(I)= [];
    end
    
    
    
    %---------------------------------------------------------------------------
    %     figure
    %     set(gcf,'position',[680   264   901   714])
    %     pcshow(colobj_mrgs{n},'Markersize',MS);
    %     hold on
    %     plot3(tdat_mrg(:,1),tdat_mrg(:,2),tdat_mrg(:,3),'.y','markersize',36)
    %     if nargin == 3
    %         plot3(eclks_mrg(:,1),eclks_mrg(:,2),eclks_mrg(:,3),'.m','markersize',32)
    %     end
    %     axis equal
    %     grid on
    %     box on
    %     FS = 12;
    %     xlabel('X','fontsize',FS,'FontName','times')
    %     ylabel('Y','fontsize',FS,'FontName','times')
    %     zlabel('Z','fontsize',FS,'FontName','times')
    %     set(gca,'FontSize',FS,'FontName','times')
    %     view([20 30])
    % saveas(gcf,['figs/ply_merged_data_attempt_v2_mrgstep',num2str(n)],'png')
        
end
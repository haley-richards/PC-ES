%-------------------------------------------------------------------------------
%
% Merge scan by merging labeled electrodes by from most intersections to
% least.
%
%-------------------------------------------------------------------------------
function [colobj_mrg,tdat_mrg] = merge_labeled_scans_tdatdirect_v2(allscans,dbg_flg)

%-------------------------------------------------------------------------------
Nsc = length(allscans);
MS  = 200;
Nel = 257;

%-------------------------------------------------------------------------------
for n = 1:Nsc
    tdats{n}   = allscans(n).tdat(1:Nel,:);
    colobjs{n} = allscans(n).colobj;
end

%-------------------------------------------------------------------------------
% Sequentially combine the point clouds
nt       = size(tdats{1},1);
tdat_mrg = NaN*ones(nt,3);
ndat_mrg =    zeros(nt,1);
nconxs   = zeros(Nsc-1,1);
id_scans = 1:Nsc;
for n = 1:(Nsc-1)
    %---------------------------------------------------------------------------
    if n == 1
        %-----------------------------------------------------------------------
        % In the initial step combine the scans with the largest number of
        % intersecting targets/stickers
        %-----------------------------------------------------------------------
        % construct the connection matrix
        cnxs = zeros(Nsc);
        for n1 = 1:(Nsc-1)
            for n2 = (n1+1):Nsc
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
        [I J]
        %-----------------------------------------------------------------------
        A   = tdats{I};
        B   = tdats{J};
    else
        %-----------------------------------------------------------------------
        % In subsequent iterations, find the largest number of intersecting
        % targets/stickers with the merged scan
        %-----------------------------------------------------------------------
        % construct the connection matrix
        cnxs = zeros(Nsc,1);
        is2  = find( isnan(tdat_mrg(:,1)) == 0);
        for n1 = 1:length(tdats)
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
    [T, rmserr]   = transform(A, B, 1);
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
    % Merge target points and point clounds
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
    % Merge point clouds!
    if n == 1
        locs_out    = (T*([colobjs{J}.Location ones(size(colobjs{J}.Location,1),1)]'))';
        locs_tmp{1} = colobjs{I}.Location;
        locs_tmp{2} = locs_out(:,1:3);
        col_tmp{1}  = colobjs{I}.Color;
        col_tmp{2}  = colobjs{J}.Color;
    else
        locs_out        = (T*([colobjs{I}.Location ones(size(colobjs{I}.Location,1),1)]'))';
        locs_tmp{end+1} = locs_out(:,1:3);
        col_tmp{end+1}  = colobjs{I}.Color;
    end
    % colobj_mrgs{n} = pointCloud(locs_mrg, 'Color', col_mrg);
    
    %---------------------------------------------------------------------------
    % Remove the individual data that was merged
    
    if n == 1
        id_mrgs{n}     = [I J];
        colobjs([I J]) = [];
        tdats([I J])   = [];
        id_scans([I J])= [];
    else
        id_mrgs{n} = id_scans(I);
        colobjs(I) = [];
        tdats(I)   = [];
        id_scans(I)= [];
    end
    
    %---------------------------------------------------------------------------
    disp(' ')
    disp('*********')
    disp(['Number of electrodes labeled and merged: ',num2str(length(find( isnan(tdat_mrg(:,1)) == 0)))])
end
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Construct a spherical mesh described by angles
naz    = 30;
nel    = 15;
angmsh = construct_sphang_msh(naz,nel,1);

%-------------------------------------------------------------------------------
% For each element of the angular mesh find the number of point cloud points
% within it and choice the scan with the most points as the one we keep
tic
for n2 = 1:length(locs_tmp)   
    rxys = sqrt(sum(double(locs_tmp{n2}(:,1:2)).^2,2));
    pts_angs{n2} = [atan2(double(locs_tmp{n2}(:,2)),double(locs_tmp{n2}(:,1))) ...
        atan(double(locs_tmp{n2}(:,3))./rxys)];
end
% Find points in each sector of each scan
np_psec = zeros(size(angmsh.angbds,1),length(locs_tmp));
locs_mrg = [];
cols_mrg = [];
for n1 = 1:size(angmsh.angbds,1)
    for n2 = 1:length(locs_tmp)   
        iss{n2} = find_points_insect(pts_angs{n2},n1,angmsh);
        np_psec(n1,n2) = length(iss{n2});
    end
    [tmp,i] = max(np_psec(n1,:));
    locs_mrg = [locs_mrg; locs_tmp{i}(iss{i},:)];
    cols_mrg = [cols_mrg; col_tmp{i}(iss{i},:)];

end
toc


MS = 100;
figure
hold on
colobj_mrg = pointCloud(locs_mrg, 'Color', cols_mrg);
pcshow(colobj_mrg,'Markersize',MS);
fac = 1.01;
plot3(fac*tdat_mrg(:,1),fac*tdat_mrg(:,2),fac*tdat_mrg(:,3),'.r','markersize',16)
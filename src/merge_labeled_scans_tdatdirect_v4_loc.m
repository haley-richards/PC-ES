%--------------------------------------------------------------------------
%
% Merge scan by merging labeled electrodes by from most intersections to
% least.
% 1. Sequentially combine the point clouds 
%
%--------------------------------------------------------------------------
function [colobj_mrg,tdat_mrg,tri_mrg,nconxs,rmserrs,omitted_scans] = merge_labeled_scans_tdatdirect_v4_loc(allscans,dbg_flg)

%--------------------------------------------------------------------------
Nsc = length(allscans);
Nel = 257;

%--------------------------------------------------------------------------
% Get cell arrays of scans and labeled electrodes 
for n = 1:Nsc
    tdats{n}   = allscans(n).tdat(1:Nel,:);
    colobjs{n} = allscans(n).colobj;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 1. Sequentially combine the point clouds 
nt       = size(tdats{1},1);        % Number of electrodes
tdat_mrg = NaN*ones(nt,3);          % Initialize the merged poins
ndat_mrg =    zeros(nt,1);          % Initialize the num of merged scans
nconxs   = zeros(Nsc*(Nsc-1)/2,1);          % Initilaize the num of connections between merged scans
id_scans = 1:Nsc;
omitted_scans = 0;
for n = 1:(Nsc-1)
    %----------------------------------------------------------------------
    if n == 1
        %------------------------------------------------------------------
        % In the initial step combine the scans with the largest number of
        % intersecting targets/stickers
        %------------------------------------------------------------------
        % construct the connection matrix
        cnxs = zeros(Nsc);
        iq   = 1;
        for n1 = 1:(Nsc-1)
            for n2 = (n1+1):Nsc
                is1 = find( isnan(tdats{n1}(:,1)) == 0);
                is2 = find( isnan(tdats{n2}(:,1)) == 0);
                cnxs(n1,n2) = length(intersect(is1,is2));
                nconxs(iq)  = cnxs(n1,n2);
                iq          = iq+1;
            end
        end

        %------------------------------------------------------------------
        % Get the location of the maximum combination
        is_max = find( max(cnxs(:)) == cnxs(:));
        [I,J]  = ind2sub(size(cnxs),is_max(1));
        % nconxs(n) = cnxs(I,J);

        %------------------------------------------------------------------
        A   = tdats{I};
        B   = tdats{J};
    else
        %------------------------------------------------------------------
        % In subsequent iterations, find the largest number of 
        % intersecting targets/stickers with the merged scan
        %------------------------------------------------------------------
        % construct the connection matrix
        cnxs = zeros(Nsc,1);
        is2  = find( isnan(tdat_mrg(:,1)) == 0);
        for n1 = 1:length(tdats)
            is1      = find( isnan(tdats{n1}(:,1)) == 0);
            cnxs(n1) = length(intersect(is1,is2));
        end
        [tmp,I] = max(cnxs);
        nconxs(n) = tmp;
        %------------------------------------------------------------------
        A   = tdat_mrg;     % A is always the prior merged scan
        B   = tdats{I};
    end

    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------

    if (rmserr < 0.015) || (n == 1)
        %----------------------------------------------------------------------
        % Merge target points and point clounds
        %   - Essentially all we are doing is averaging points together
        %   - we make sure we don't try to average a number with a NaN
        %   - we make sure that each number is weighted equally, i.e. if in the
        %     merged scan there was already averaging of a given point done, we
        %     account for this and don't just weight it equally with the next
        %     point.
        %----------------------------------------------------------------------
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
            tri_tmp{1}  = allscans(I).tri;
            tri_tmp{2}  = allscans(J).tri;
        else
            locs_out        = (T*([colobjs{I}.Location ones(size(colobjs{I}.Location,1),1)]'))';
            locs_tmp{end+1} = locs_out(:,1:3);
            col_tmp{end+1}  = colobjs{I}.Color;
            tri_tmp{end+1}  = allscans(I).tri;
        end
        % colobj_mrgs{n} = pointCloud(locs_mrg, 'Color', col_mrg);
        
        %----------------------------------------------------------------------
        if dbg_flg == 1
            disp(' ')
            disp('*********')
            disp(['Number of electrodes labeled and merged: ',num2str(length(find( isnan(tdat_mrg(:,1)) == 0)))])
        end
    else
        if n == 1
            error('Case not considered, presumptively 1st step has an outlier')
        end
        disp('Scan fit too bad: Not adding, removing the preposed scan fit')
        omitted_scans = omitted_scans+1;
    end
    %------------------------------------------------------------------
    % Remove the individual data that was merged
    if n == 1
        id_mrgs{n}     = [I J];
        colobjs([I J]) = [];
        allscans([I J])= [];
        tdats([I J])   = [];
        id_scans([I J])= [];
    else
        id_mrgs{n} = id_scans(I);
        colobjs(I) = [];
        allscans(I)= [];
        tdats(I)   = [];
        id_scans(I)= [];
    end

end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Construct a spherical mesh described by angles
naz    = 30;
nel    = 15;
angmsh = construct_sphang_msh(naz,nel,dbg_flg);

%--------------------------------------------------------------------------
% For each element of the angular mesh find the number of point cloud points
% within it and choice the scan with the most points as the one we keep
% tic
for n2 = 1:length(locs_tmp)
    rxys = sqrt(sum(double(locs_tmp{n2}(:,1:2)).^2,2));
    pts_angs{n2} = [atan2(double(locs_tmp{n2}(:,2)),double(locs_tmp{n2}(:,1))) ...
        atan(double(locs_tmp{n2}(:,3))./rxys)];
end
%-------------------------------------------------------------------------------
% Find points in each sector of each scan. Collect all the indices from
% each scan that will be kept before doing any combining of scans
np_psec = zeros(size(angmsh.angbds,1),length(locs_tmp));
locs_mrg = [];
cols_mrg = [];
istot    = cell(1,length(locs_tmp));
for n1 = 1:size(angmsh.angbds,1)
    for n2 = 1:length(locs_tmp)
        iss{n2} = find_points_insect(pts_angs{n2},n1,angmsh);
        np_psec(n1,n2) = length(iss{n2});
    end
    [tmp,i] = max(np_psec(n1,:));
    istot{i} = [istot{i}; iss{i}];
end
%--------------------------------------------------------------------------
% After we got all the indices we want to keep from each scan, then we
% combine point clouds, and triangulations
if dbg_flg == 1
    for i = 1:length(locs_tmp)
        disp([num2str(i),': p and t make sense: ',num2str([size(locs_tmp{i},1) max(tri_tmp{i}(:))])])
    end
end
tri_mrg = [];
for i = 1:length(locs_tmp)
    if ~isempty(istot{i})
        ris   = setdiff((1:size(locs_tmp{i},1))',istot{i});
        [p,t] = update_tris_due_to_rmvnodes(locs_tmp{i},tri_tmp{i},ris);
        tri_mrg  = [tri_mrg; t + size(locs_mrg,1)];
        col_tmp{i}(ris,:) = [];
        locs_mrg = [locs_mrg; p];
        cols_mrg = [cols_mrg; col_tmp{i}];
    end
end
% toc
colobj_mrg = pointCloud(locs_mrg, 'Color', cols_mrg);

if dbg_flg == 1
    MS = 100;
    figure
    hold on
    pcshow(colobj_mrg,'Markersize',MS);
    fac = 1.01;
    plot3(fac*tdat_mrg(:,1),fac*tdat_mrg(:,2),fac*tdat_mrg(:,3),'.r','markersize',16)
    figure
    trisurf(tri_mrg,locs_mrg(:,1),locs_mrg(:,2),locs_mrg(:,3),'facecolor','cyan','linestyle','none')
    camlight left
end
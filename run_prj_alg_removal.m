%-------------------------------------------------------------------------------
%
% Run the projection algorithm
%
%-------------------------------------------------------------------------------
function [is1c, pc_out] = run_prj_alg_removal(colobj,pts0)
% clear
% clc
% close all
% load testdat


%-------------------------------------------------------------------------------
% Take the points 3D and only keep those above the z = 0 mark
locs0 = double(colobj.Location);
cols0 = colobj.Color;
is1   = find( locs0(:,3) > 0);
locs  = locs0(is1,:);
cols  = cols0(is1,:);
if size(pts0,1) > 0
    is1b   = find( pts0(:,3) > 0);
    pts   = pts0(is1b,:);
end

%-------------------------------------------------------------------------------
% Take 90% of the convex hull to do image processing analysis on
K    = convhull(locs(:,1),locs(:,2));
in   = inpolygon(locs(:,1),locs(:,2),0.9*locs(K,1),0.9*locs(K,2));
is2   = find(in == 1);
locs = locs(is2,:);
cols = cols(is2,:);
cols = mean(cols,2);

%-------------------------------------------------------------------------------
% Construct a meshgrid of the locations values
npx = 1000;
xs      = linspace(min(locs(:,1)),max(locs(:,1)),npx);
ys      = linspace(min(locs(:,2)),max(locs(:,2)),npx);
[xm,ym] = meshgrid(xs,ys);
xys     = [xm(:) ym(:)];

%-------------------------------------------------------------------------------
% Interpolate onto the meshgrid
F     = scatteredInterpolant(locs(:,1:2),cols,'linear','none');
img_v = F(xys);
img_m = reshape(img_v,size(xm,1),size(xm,2));
%-------------------------------------------------------------------------------
% Interpolate the points into the image frame
if size(pts0,1) > 0
    pts_m(:,1) = interp1(xs,1:npx,pts(:,1));
    pts_m(:,2) = interp1(ys,1:npx,pts(:,2));
else
    pts_m = [];
end

%     imagesc(img_m)
%     set(gca, 'YDir','normal')
%     colormap(gray)
%     axis equal
%     drawnow

%-------------------------------------------------------------------------------
% Run circle center finding algorithm
[ris, cents] = find_circs_to_remove(img_m,pts_m);
if size(cents, 1) > 0
    %---------------------------------------------------------------------------
    % convert the centers from pixel space to projected 2D space
    cent_m      = zeros(size(cents,1),2);
    cent_m(:,1) = interp1(1:npx,xs,cents(:,1));
    cent_m(:,2) = interp1(1:npx,ys,cents(:,2));
    
    %-------------------------------------------------------------------------------
    % Transform the centers to 3D space, i.e. find
    newtdats = zeros(size(cents,1),3);
    n_colpts=80;
    locs=colobj.Location;
    cols_out=colobj.Color;
    for n = 1:size(cents,1)
        %---------------------------------------------------------------------------
        [tmp,j] = min( (locs(:,1) - cent_m(n,1)).^2 + (locs(:,2) - cent_m(n,2)).^2);
         [locs_sort, idx]=sort(((locs(:,1) - cent_m(n,1)).^2 + (locs(:,2) - cent_m(n,2)).^2), 'ascend');
        cols_out(idx(1:n_colpts), :)=repmat([0, 0, 0], n_colpts, 1); 
    end
end
pc_out=pointCloud(locs, 'Color', cols_out);
is1c  = is1b(ris);


%-------------------------------------------------------------------------------
%
% Run circle center finding algorithm
%
%-------------------------------------------------------------------------------
function [ris, cents] = find_circs_to_remove(img_m,pts_m,dbg_flg)
% save testdata
% clear

% clear
% clc
% close all
% load testdata

%---------------------------------------------------------------------------
% Plot the image and pick an electrode as the template
figure
set(gcf,'position',[257         124        1198         825])
imagesc(img_m)
hold on
if size(pts_m,1) > 0
    plot(pts_m(:,1),pts_m(:,2),'.r','markersize',12)
end
colormap gray
hold on
axis equal
set(gca, 'YDir','normal')
title('Pick incorrectly marked electrode and Right click to finish')
butt    = 1;
ris = [];
cent_count=1;
while (butt == 1) 
    [x1,y1,butt] = ginputWhite(1);
    cents(cent_count, :)=[x1, y1];
    if butt == 1
        [tmp,i]      = min( sum( (pts_m - repmat([x1 y1],size(pts_m,1),1)).^2,2));
        cent_count=cent_count+1;
        ris          = [ris; i];
        pts_m(ris,:) = NaN;
        clf
        imagesc(img_m)
        hold on
        if size(pts_m,1) > 0
            plot(pts_m(:,1),pts_m(:,2),'.r','markersize',18)
        end               
        colormap gray
        hold on
        axis equal
        set(gca, 'YDir','normal')        
    end
end






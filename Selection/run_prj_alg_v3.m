%-------------------------------------------------------------------------------
%
% Run the projection algorithm
%
%-------------------------------------------------------------------------------
function [newtdats] = run_prj_alg_v3(colobj,pts0,ttl_str,dbg_flg)
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
% cols = mean(cols,2);

%-------------------------------------------------------------------------------
% Construct a meshgrid of the locations values
npx = 1000;
xs      = linspace(min(locs(:,1)),max(locs(:,1)),npx);
ys      = linspace(min(locs(:,2)),max(locs(:,2)),npx);
[xm,ym] = meshgrid(xs,ys);
xys     = [xm(:) ym(:)];

%-------------------------------------------------------------------------------
% Interpolate onto the meshgrid
img_m = zeros(size(xm,1),size(xm,2),3);
for i = 1:3
    F     = scatteredInterpolant(locs(:,1:2),double(cols(:,i))/255,'linear','none');
    img_v = F(xys);
    img_m(:,:,i) = reshape(img_v,size(xm,1),size(xm,2));
end
%-------------------------------------------------------------------------------
% Interpolate the points into the image frame
if size(pts0,1) > 0
    pts_m(:,1) = interp1(xs,1:npx,pts(:,1));
    pts_m(:,2) = interp1(ys,1:npx,pts(:,2));
else
    pts_m = [];
end
if dbg_flg == 1
    imagesc(img_m)
    set(gca, 'YDir','normal')
    % colormap(gray)
    axis equal
    drawnow
end
%-------------------------------------------------------------------------------
% Run circle center finding algorithm
[cents] = find_circ_cents_v3(img_m,pts_m,ttl_str,dbg_flg);

if size(cents,1) > 0
    %---------------------------------------------------------------------------
    % convert the centers from pixel space to projected 2D space
    cent_m      = zeros(size(cents,1),2);
    cent_m(:,1) = interp1(1:npx,xs,cents(:,1));
    cent_m(:,2) = interp1(1:npx,ys,cents(:,2));
    
    %-------------------------------------------------------------------------------
    % Transform the centers to 3D space, i.e. find
    newtdats = zeros(size(cents,1),3);
    for n = 1:size(cents,1)
        %---------------------------------------------------------------------------
        [tmp,j] = min( (locs(:,1) - cent_m(n,1)).^2 + (locs(:,2) - cent_m(n,2)).^2);
        newtdats(n,:) = locs0(is1(is2(j)),:);
    end
    
else
    newtdats = [];
end
% axes(handles.mainax);
% hold on
% plot3(newtdats(:,1),newtdats(:,2),newtdats(:,3),'.g','markersize',24)





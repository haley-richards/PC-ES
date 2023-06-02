function [newtdats, colobj_out] = run_prj_alg_v2_sub(colobj,pts0,dbg_flg, colobj0)
% clear
% clc
% close all
% load testdat

locs_out=colobj0.Location;
cols_out=colobj0.Color;
%-------------------------------------------------------------------------------
% Take the points 3D and only keep those above the z = 0 mark
locs0 = double(colobj.Location);
cols0 = colobj.Color;
is1   = find( locs0(:,3) > 0);
locs  = locs0;
cols  = cols0;
if size(pts0,1) > 0
    is1b   = find( pts0(:,3) > 0);
    pts   = pts0(is1b,:);
end

%-------------------------------------------------------------------------------
% Take 90% of the convex hull to do image processing analysis on
K    = convhull(locs(:,1),locs(:,2));
in   = inpolygon(locs(:,1),locs(:,2),0.9*locs(K,1),0.9*locs(K,2));
is2   = find(in == 1);
locs = locs;
cols = cols;
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
 img_m1 = lin2rgb(img_m);
  
   % chart = colorChecker(img_m1);
%     displayChart(chart)
%     registrationPoints = chart_sRGB.RegistrationPoints;
    
    %illuminant_groundtruth = measureIlluminant(chart)
    img_m2=img_m;
    img_m2(isnan(img_m2))=0;
    img_m2=double(img_m2);
    percentileToExclude = 0;
illuminant_wp1 = illumwhite(img_m2,percentileToExclude);
img_m3 = chromadapt(img_m2,illuminant_wp1,'ColorSpace','linear-rgb');
  imagesc(img_m3)
    set(gca, 'YDir','normal')
    % colormap(gray)
    axis equal
    drawnow

%-------------------------------------------------------------------------------
% Run circle center finding algorithm
[cents] = find_circ_cents_hr(img_m3,pts_m,dbg_flg);

%-------------------------
%change area around point to be red
% tdat_sub=pc_in.Location(tdat_idx, :);
% new_cols=pc_in.Color;
%colobj.Color(pc_sub(idx(1:10)), :)=repmat([250, 0, 0], 10, 1);
if size(cents, 1) > 0
    %---------------------------------------------------------------------------
    % convert the centers from pixel space to projected 2D space
    cent_m      = zeros(size(cents,1),2);
    cent_m(:,1) = interp1(1:npx,xs,cents(:,1));
    cent_m(:,2) = interp1(1:npx,ys,cents(:,2));
    
    %-------------------------------------------------------------------------------
    % Transform the centers to 3D space, i.e. find
    newtdats = zeros(size(cents,1),3);
    n_colpts=50;
    for n = 1:size(cents,1)
        %---------------------------------------------------------------------------
        [tmp,j] = min( (locs(:,1) - cent_m(n,1)).^2 + (locs(:,2) - cent_m(n,2)).^2);
         [locs_sort, idx]=sort(((locs(:,1) - cent_m(n,1)).^2 + (locs(:,2) - cent_m(n,2)).^2), 'ascend');
        cols_out(idx(1:n_colpts), :)=repmat([250, 0, 0], n_colpts, 1); 
        newtdats(n,:) = locs0((j),:);
    end
    
else
    newtdats = [];
end
colobj_out=pointCloud(locs_out, 'Color', cols_out);
close
end




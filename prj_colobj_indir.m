%-------------------------------------------------------------------------------
%     
% Project onto the nth direction and produce the image
% 
%-------------------------------------------------------------------------------
function prj_colobj_indir(colobj,tdat,i,angmsh,dbg_flg)

%-------------------------------------------------------------------------------
% To be consistent with recent work we want to rotate coordinate frame so
% the direction of interest align with the +z-axis. So we can do a z-axis
% rotation back to the +y-axis and then a y-axis rotation up the +z-axis
az = angmsh.angcnts(i,1);
el = angmsh.angcnts(i,2);
Rz = Rotz(az);
Ry = Roty(pi/2-el);
Rz2= Rotz(pi/2);

%-------------------------------------------------------------------------------
% perform the rotation
locs = double(colobj.Location);
cols = double(colobj.Color);
% cols = mean(cols,2);
locs = (Rz2'*Ry'*Rz'*(locs'))';
tdat = (Rz2'*Ry'*Rz'*(tdat'))';


% figure
% colobj = pointCloud(locs(:,1:3), 'Color', colobj.Color);
% pcshow(colobj)
% 
% clear
%-------------------------------------------------------------------------------
% Only keep positive z-valued points
is   = find( locs(:,3) > 0);
locs = locs(is,:);
cols = cols(is,:);

    
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
for n = 1:3
    F     = scatteredInterpolant(locs(:,1:2),cols(:,n)/255,'linear','none');
    img_v = F(xys);
    img_m(:,:,n) = reshape(img_v,size(xm,1),size(xm,2));
end
%-------------------------------------------------------------------------------
% Interpolate the points into the image frame
if size(tdat,1) > 0
    pts_m(:,1) = interp1(xs,1:npx,tdat(:,1));
    pts_m(:,2) = interp1(ys,1:npx,tdat(:,2));
else
    pts_m = [];
end
if dbg_flg == 1
    figure
    imagesc(img_m)    
    set(gca, 'YDir','normal')
    colormap(gray)
    axis equal
    drawnow
end

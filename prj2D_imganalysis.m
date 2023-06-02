function pts3D = prj2D_imganalysis(colobj,dbg_flg)

%-------------------------------------------------------------------------------
% Interpolate onto a 2D grid: assume the location matrix is in meters and
% convert them to 1/2^2 mm's
xlims   = double(colobj.XLimits*1000);
ylims   = double(colobj.YLimits*1000);
xs      = linspace(xlims(1),xlims(2),4*round(xlims(2)-xlims(1)));
ys      = linspace(ylims(1),ylims(2),4*round(ylims(2)-ylims(1)));
[xm,ym] = meshgrid(xs,ys);
F       = scatteredInterpolant(double(colobj.Location(:,1)*1000),double(colobj.Location(:,2)*1000), ...
    mean(double(colobj.Color),2),'linear','linear');
xyps    = [xm(:) ym(:)];
Vq      = F(xyps);
%-------------------------------------------------------------------------------
% Scale from 0 to 1 and interpolate on an image grid
isneg     = find(Vq < 0);
Vq(isneg) = 0;
img2d     = reshape(Vq/max(Vq),size(xm,1),size(xm,2));

%-------------------------------------------------------------------------------
% if dbg_flg == 1
%     figure
%     imagesc(img2d)
%     colormap gray
% end

%-------------------------------------------------------------------------------
% Find circles
[centers,radii] = imfindcircles(img2d,[5 12],'ObjectPolarity','bright');% ,'Sensitivity',0.95);
% [centersBright,radiiBright,metricBright] = imfindcircles(img2d,[6 12], ...
%     'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.1);

if dbg_flg == 1
    figure
    imagesc(img2d)
    colormap gray
    hold on
    h = viscircles(centers,radii);
    axis equal
end

%-------------------------------------------------------------------------------
% Convert the 2D centers to 3D points
centers = round(centers);
pts     = [xs(centers(:,1))' ys(centers(:,2))'];
locs    = double(colobj.Location*1000);
pts3D   = zeros(size(pts,1),3);
for n = 1:size(pts,1)
    ds         = (locs(:,1)-pts(n,1)).^2 + (locs(:,2)-pts(n,2)).^2;
    [tmp,i]    = min(ds);
    pts3D(n,:) = locs(i,:); 
end
pts3D = pts3D/1000;

return


img2d_bw     = imbinarize(img2d,'adaptive','Sensitivity',0.75);
% % imshow(bw)
% figure
% imagesc(img2d_bw)
% colormap gray
% 
% % figure
% % img2d_bwb = imfill(img2d_bw,'holes');
% % imshow(img2d_bwb)
% 
% [B,L] = bwboundaries(img2d_bw,'noholes');
% imshow(label2rgb(L,@jet,[.5 .5 .5]))
% hold on
% for k = 1:length(B)
%   boundary = B{k};
%   plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
% end


[centersBright,radiiBright,metricBright] = imfindcircles(img2d_bw,[6 12], ...
    'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.1);
figure
imagesc(img2d_bw)
colormap gray
hold on
h = viscircles(centers,radii);

return
% d = drawline;
% pos = d.Position;
% diffPos = diff(pos);
% diameter = hypot(diffPos(1),diffPos(2))

[centers,radii] = imfindcircles(img2d_bw,[5 12],'ObjectPolarity','bright');% ,'Sensitivity',0.95);
% [centers,radii] = imfindcircles(img2d,[5 12],'ObjectPolarity','bright','Method','twostage');
hold on
h = viscircles(centers,radii);



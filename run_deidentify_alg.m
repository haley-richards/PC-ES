%--------------------------------------------------------------------------
%
% Run the projection algorithm
%
%--------------------------------------------------------------------------
function handles = run_deidentify_alg(handles)


%--------------------------------------------------------------------------
colobj = update_PLY_pts(handles);

%--------------------------------------------------------------------------
% Get new points
dbg_flg = 0;
di_is   = run_prj_deident_alg(colobj,dbg_flg);
%--------------------------------------------------------------------------
% Make these points gray
cols           = handles.colobj.Color;
cols(di_is,:)  = repmat([200 200 200],length(di_is),1);
colobj0        = pointCloud(handles.colobj.Location, 'Color', cols);
handles.colobj = colobj0;

%--------------------------------------------------------------------------
figure
set(gcf,'position',[258         213        1025         739])
hold on
az2        = str2double(get(handles.cur_az,'String'))*pi/180;
el2        = str2double(get(handles.cur_el,'String'))*pi/180;
[colobj,tdat,eclks,cent0] = update_PLY_pts(handles,[az2 el2]);
MS = 100;
pcshow(colobj,'Markersize',MS);
axis equal
view(2)




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% 
%
%--------------------------------------------------------------------------
function di_is = run_prj_deident_alg(colobj,dbg_flg)


%--------------------------------------------------------------------------
% Take the points 3D and only keep those above the z = 0 mark
locs0 = double(colobj.Location);
cols0 = colobj.Color;
is1   = find( locs0(:,3) > 0);
locs  = locs0(is1,:);
cols  = cols0(is1,:);

%--------------------------------------------------------------------------
% Construct a meshgrid of the locations values
npx = 1000;
xs      = linspace(min(locs(:,1)),max(locs(:,1)),npx);
ys      = linspace(min(locs(:,2)),max(locs(:,2)),npx);
[xm,ym] = meshgrid(xs,ys);
xys     = [xm(:) ym(:)];

%--------------------------------------------------------------------------
% Interpolate onto the meshgrid
img_m = zeros(size(xm,1),size(xm,2),3);
for i = 1:3
    F     = scatteredInterpolant(locs(:,1:2),double(cols(:,i))/255,'linear','none');
    img_v = F(xys);
    img_m(:,:,i) = reshape(img_v,size(xm,1),size(xm,2));
end

%--------------------------------------------------------------------------
% Run circle center finding algorithm
[cents] = find_boundary_identifying_region(img_m);

if size(cents,1) > 0
    %----------------------------------------------------------------------
    % convert the centers from pixel space to projected 2D space
    cent_m      = zeros(size(cents,1),2);
    cent_m(:,1) = interp1(1:npx,xs,cents(:,1));
    cent_m(:,2) = interp1(1:npx,ys,cents(:,2));
    
    %----------------------------------------------------------------------
    % Find points within the boundary selected    
    IN    = inpolygon(locs(:,1), locs(:,2), cent_m(:,1), cent_m(:,2));
    di_is = is1(find( IN == 1));        
    
else
    di_is = [];
end
% axes(handles.mainax);
% hold on
% plot3(newtdats(:,1),newtdats(:,2),newtdats(:,3),'.g','markersize',24)


%--------------------------------------------------------------------------
%
% Run circle center finding algorithm
%
%--------------------------------------------------------------------------
function [centers] = find_boundary_identifying_region(img_m)

%---------------------------------------------------------------------------
% Plot the image and pick an electrode as the template
figure
set(gcf,'position',[257         124        1198         825])
imagesc(img_m)
colormap gray
hold on
axis equal
set(gca, 'YDir','normal')
title('Pick in a counter clockwise direction the identifying region of the scan, 1 = add, 2 = fix, 3 = done')
butt    = 1;
centers = [];
while (butt == 1) || (butt == 2)
    [x1,y1,butt] = ginputWhite(1);
    if butt == 1        
        plot(x1,y1,'.g','markersize',18)
        centers = [centers; x1,y1];
    elseif butt == 2
        [tmp,i] = min( sum( (centers - repmat([x1 y1],size(centers,1),1)).^2,2));
        centers(i,:) = [];
        clf
        imagesc(img_m)
        hold on
        if size(centers,1) > 0 
            plot(centers(:,1),centers(:,2),'.g','markersize',18)
        end        
        % colormap gray
        hold on
        axis equal
        set(gca, 'YDir','normal')        
    end
end





%-------------------------------------------------------------------------------
% 
% pts_angs are Nx2, 
%   column 1 is azimuth angles (-pi to pi or 0 to 2*pi)
%   column 2 is elevation angles (-pi/2 to pi/2)
% 
%-------------------------------------------------------------------------------
function is = find_points_insect(pts_angs,n,angmsh)

%-------------------------------------------------------------------------------
is = find( ...
    ( ((angmsh.angbds(n,1) < pts_angs(:,1)) & (pts_angs(:,1) < angmsh.angbds(n,2) )) | ...
    ((angmsh.angbds(n,1) < pts_angs(:,1)+2*pi) & (pts_angs(:,1)+2*pi < angmsh.angbds(n,2) )) | ...
    ((angmsh.angbds(n,1) < pts_angs(:,1)-2*pi) & (pts_angs(:,1)-2*pi < angmsh.angbds(n,2) ) )) & ...
    (angmsh.angbds(n,3) < pts_angs(:,2) ) & (pts_angs(:,2) < angmsh.angbds(n,4) ) );


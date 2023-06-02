%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Update the figure
function [colobj,tdat,eclks,cent0] = update_PLY_pts_notGUI(handles)

tdat  = [];
eclks = [];
%-------------------------------------------------------------------------------
% If there is an origin defined subtract off the origin
cent = [ ...
    handles.org_x ...
    handles.org_y ...
    handles.org_z];
%-------------------------------------------------------------------------------
% Load the surface and stop if there is no STL
if isfield(handles,'colobj') == 1
    colobj = handles.colobj;
    locs   = colobj.Location;
    cols   = colobj.Color;
    cent0  = mean(locs,1);
    locs   = locs - repmat(cent0,size(locs,1),1);
else
    error('stop')
end
%-------------------------------------------------------------------------------
% If there is an origin defined subtract off the origin
if isfield(handles,'tdat') == 1
    tdat = handles.tdat;
    tdat = tdat - repmat(cent0,size(tdat,1),1);
end
if isfield(handles,'eclks') == 1
    eclks = handles.eclks;
    eclks = eclks - repmat(cent0,size(eclks,1),1);
end


%-------------------------------------------------------------------------------
% Rotate the STL to the initial frame
az = handles.ini_rz*pi/180;
el = handles.ini_rx*pi/180;
Rx          = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
if isfield(handles,'tdat') == 1
    tdat        = (Rx*(tdat'))';
end
if isfield(handles,'eclks') == 1
    eclks       = (Rx*(eclks'))';
end
locs        = (Rx*(locs'))';
Rz          = [cos(az) sin(az) 0;-sin(az) cos(az) 0; 0 0 1];
if isfield(handles,'tdat') == 1
    tdat        = (Rz*(tdat'))';
end
if isfield(handles,'eclks') == 1
    eclks       = (Rz*(eclks'))';
end
locs        = (Rz*(locs'))';

%-------------------------------------------------------------------------------
locs    = locs    - repmat(cent,size(locs,1),1);
if isfield(handles,'tdat') == 1
    tdat  = tdat  - repmat(cent,size(tdat,1),1);
end
if isfield(handles,'eclks') == 1
    eclks = eclks - repmat(cent,size(eclks,1),1);
end

%---------------------------------------------------------------------------
% Rotate the STL to the temporary frame
az2 = handles.cur_az*pi/180;
el2 = handles.cur_el*pi/180;
Rz          = [cos(az2) sin(az2) 0;-sin(az2) cos(az2) 0; 0 0 1];
if isfield(handles,'tdat') == 1
    tdat        = (Rz*(tdat'))';
end
if isfield(handles,'eclks') == 1
    eclks       = (Rz*(eclks'))';
end
locs          = (Rz*(locs'))';
Rx          = [1 0 0; 0 cos(el2) -sin(el2); 0 sin(el2) cos(el2)];
if isfield(handles,'tdat') == 1
    tdat        = (Rx*(tdat'))';
end
if isfield(handles,'eclks') == 1
    eclks       = (Rx*(eclks'))';
end
locs          = (Rx*(locs'))';
 
%---------------------------------------------------------------------------
% If the plot plane flag is selected, we don't actually plot a plane,
% because its a point cloud and not a surface triangulation it messes up
% the visualization. So, instead we remove all points below the plane,
% which is just about like plotting an invisible wall..
if handles.plotplane == 1
    z0  = handles.zpl_z;
    ris = find( locs(:,3) < z0);
    locs(ris,:) = [];
    cols(ris,:) = [];
end
%---------------------------------------------------------------------------
colobj        = pointCloud(locs, 'Color', cols);    

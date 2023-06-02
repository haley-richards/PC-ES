%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Update the figure
function [colobj,tdat,eclks,cent0] = update_PLY_pts(handles,azels2)

tdat  = [];
eclks = [];
%-------------------------------------------------------------------------------
% If there is an origin defined subtract off the origin
cent = [0 0 0]; %  ...
%     str2double(get(handles.org_x,'String')) ...
%     str2double(get(handles.org_y,'String')) ...
%     str2double(get(handles.org_z,'String'))];
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
az = handles.ini_rz; % 0;% str2double(get(handles.ini_rz,'String'))*pi/180;
el = handles.ini_rx; % 90*pi/180; % str2double(get(handles.ini_rx,'String'))*pi/180;
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
if nargin < 2
    az2 = str2double(get(handles.cur_az,'String'))*pi/180;
    el2 = str2double(get(handles.cur_el,'String'))*pi/180;
else % Use the input azimum and elevation angles (assumed in radians)
    az2 = azels2(1);
    el2 = azels2(2);
end
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
    z0  = str2double(get(handles.zpl_z,'String'));
    ris = find( locs(:,3) < z0);
    locs(ris,:) = [];
    cols(ris,:) = [];
end
%---------------------------------------------------------------------------
colobj        = pointCloud(locs, 'Color', cols);    

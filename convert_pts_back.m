%-------------------------------------------------------------------------------
%
% This function uses the information in the GUI to convert the points back to
% the original frame
%
%-------------------------------------------------------------------------------
function pts0 = convert_pts_back(pts,handles,cent0,azels2)

if size(pts,1) ~= 3
    pts = pts';
end

%-------------------------------------------------------------------------------
cent = [0 0 0]'; %  ...
%     str2double(get(handles.org_x,'String')) ...
%     str2double(get(handles.org_y,'String')) ...
%     str2double(get(handles.org_z,'String'))]';
%---------------------------------------------------------------------------
xbs = [-0.1 0.1]; % [str2double(get(handles.zpl_x1,'String')) ...
    % str2double(get(handles.zpl_x2,'String'))];
ybs = [-0.1 0.1]; % [str2double(get(handles.zpl_y1,'String')) ...
    % str2double(get(handles.zpl_y2,'String'))];
z0  = 0;% str2double(get(handles.zpl_z,'String'));

%---------------------------------------------------------------------------
% Convert the point back to the initial coordinates
if nargin < 4
    az2 = str2double(get(handles.cur_az,'String'))*pi/180;
    el2 = str2double(get(handles.cur_el,'String'))*pi/180;
else % Use the input azimum and elevation angles (assumed in radians)
    az2 = azels2(1);
    el2 = azels2(2);
end
Rx   = [1 0 0; 0 cos(el2) -sin(el2); 0 sin(el2) cos(el2)];
newp = (Rx')*pts;
Rz   = [cos(az2) sin(az2) 0;-sin(az2) cos(az2) 0; 0 0 1];
newp = (Rz')*newp;
newp = newp + repmat(cent,1,size(newp,2));
%-------------------------------------------------------
% Rotate the STL to the initial frame
az   = 0;% str2double(get(handles.ini_rz,'String'))*pi/180;
el   = 90*pi/180; % str2double(get(handles.ini_rx,'String'))*pi/180;
Rz   = [cos(az) sin(az) 0;-sin(az) cos(az) 0; 0 0 1];
newp = (Rz')*newp;
Rx   = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
newp = (Rx')*newp;
%-------------------------------------------------------
newp = newp + repmat(cent0',1,size(newp,2));
pts0 = newp';


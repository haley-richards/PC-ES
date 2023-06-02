%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Update the figure
function [pin,cent0] = rotate_general_pts(pin,handles,azels2)

%-------------------------------------------------------------------------------
% If there is an origin defined subtract off the origin
cent = [0 0 0];%  ...
%     str2double(get(handles.org_x,'String')) ...
%     str2double(get(handles.org_y,'String')) ...
%     str2double(get(handles.org_z,'String'))];
%-------------------------------------------------------------------------------
% Load the surface and stop if there is no STL
if isfield(handles,'colobj') == 1
    colobj = handles.colobj;
    locs   = colobj.Location;
    cent0  = mean(locs,1);    
else
    error('stop')
end
pin = pin - repmat(cent0,size(pin,1),1);

%-------------------------------------------------------------------------------
% Rotate the STL to the initial frame
az = 0; % str2double(get(handles.ini_rz,'String'))*pi/180;
el = 90*pi/180; % str2double(get(handles.ini_rx,'String'))*pi/180;
Rx          = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
pin         = (Rx*(pin'))';
Rz          = [cos(az) sin(az) 0;-sin(az) cos(az) 0; 0 0 1];
pin         = (Rz*(pin'))';

%-------------------------------------------------------------------------------
pin    = pin    - repmat(cent,size(pin,1),1);

%---------------------------------------------------------------------------
% Rotate the STL to the temporary frame
if nargin < 3
    az2 = str2double(get(handles.cur_az,'String'))*pi/180;
    el2 = str2double(get(handles.cur_el,'String'))*pi/180;
else % Use the input azimum and elevation angles (assumed in radians)
    az2 = azels2(1);
    el2 = azels2(2);
end
Rz          = [cos(az2) sin(az2) 0;-sin(az2) cos(az2) 0; 0 0 1];
pin         = (Rz*(pin'))';
Rx          = [1 0 0; 0 cos(el2) -sin(el2); 0 sin(el2) cos(el2)];
pin         = (Rx*(pin'))';
 

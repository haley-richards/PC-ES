
function varargout = PCES(varargin)
% PCES MATLAB code for PCES.fig
%      PCES, by itself, creates a new PCES or raises the existing
%      singleton*.
%
%      H = PCES returns the handle to a new PCESR or the handle to
%      the existing singleton*.
%
%      PCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCES.M with the given input arguments.
%
%      PCES('Property','Value',...) creates a new PCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCES_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCES_Fcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCES

% Last Modified by GUIDE v2.5 02-Oct-2022 00:12:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PCES_OpeningFcn, ...
    'gui_OutputFcn',  @PCES_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PCES is made visible.
function PCES_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCES (see VARARGIN)

% Choose default command line output for PCES
hObject.Name='PCES';
handles.output = hObject;
%add correct paths and home directory
addpath(genpath(cd))
%load homepage
im=imread('homepage2.jpg');
imshow(im)
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = PCES_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figuref
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pick_vx_stl.
function pick_vx_stl_Callback(hObject, eventdata, handles)
% hObject    handle to pick_vx_stl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath(genpath('S:/digihisto/Ethan/EKM_utility'))
% addpath mfiles
% addpath(genpath('/Volumes/jumbo/digihisto/Ethan/EKM_utility'))
% addpath mfiles
% addpath(genpath('/jumbo/digihisto/Ethan/EKM_utility'))


%-------------------------------------------------------------------------------
% An STL file to process
[file,path] = uigetfile('*.PLY');

%-------------------------------------------------------------------------------
colobj = pcread([path,'/',file]);
set(handles.file_name_txtbx,'String',file)
%-------------------------------------------------------------------------------
% Remove the black points
ris          = find( sum(double(colobj.Color),2) < 2);
locs         = colobj.Location;
locs(ris,:)  = [];
cols         = colobj.Color;
cols(ris,:)  = [];
colobj       = pointCloud(locs, 'Color', cols);
%-------------------------------------------------------------------------------
% Choose default command line output for viewMRdat
handles.colobj    = colobj;
handles.file      = file;
handles.cropped   = 0;
handles.plotplane = 0;
handles.ini_rz    = 0;         
handles.ini_rx    = 90*pi/180; 
handles.org_x     = 0;       
handles.org_y     = 0;       
handles.org_z     = 0;
handles.zpl_x1    = -0.1;
handles.zpl_x2    =  0.1;
handles.zpl_y1    = -0.1;
handles.zpl_y2    =  0.1;
handles.zpl_z     = 0;
update_ptcld_plot(hObject,handles)
%-------------------------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in crop_ply_dat.
function crop_ply_dat_Callback(hObject, eventdata, handles)
% hObject    handle to crop_ply_dat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colobj = handles.colobj;
handles.plotplane = 0;  % Plot plane cannot be selected

%--------------------------------------------------------------------------
% We want to go through 3 views of the PLY data and crop all the needed
% data
axes(handles.mainax);
for n = 1:2
    %----------------------------------------------------------------------
    locs = colobj.Location;
    cols = colobj.Color;
    %----------------------------------------------------------------------
    % Rotate the points
    if n == 1
        R = eye(3);
    elseif n == 2
        t = 90*pi/180;
        R = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
    end
    locs = (R*(locs'))';
    colobj         = pointCloud(locs, 'Color', cols);
    handles.colobj = colobj;
    update_ptcld_plot(hObject,handles)
    
    %----------------------------------------------------------------------
    view(2)
    set(handles.message_texts,'String','Pick the Top Left Extent of Interest')
    [xleft,ytop] = ginputWhite(1);
    view(2)
    set(handles.message_texts,'String','Pick the Bottom Right Extent of Interest')
    [xright,ybot] = ginputWhite(1);
    ris    = find( (colobj.Location(:,1) < xleft) | (colobj.Location(:,1) > xright) | ...
        (colobj.Location(:,2) < ybot) | (colobj.Location(:,2) > ytop) );
    locs(ris,:) = [];
    cols(ris,:) = [];
    colobj      = pointCloud(locs, 'Color', cols);
    handles.colobj = colobj;
    update_ptcld_plot(hObject,handles)
    pause(1)
    %---------------------------------------------------------------------------
    % Un-Rotate the points
    locs = ((R')*(locs'))';
    colobj         = pointCloud(locs, 'Color', cols);
    handles.colobj = colobj;
            
end
handles.cropped = 1;
handles.colobj0   = colobj;
handles.tdat=NaN(256, 3);
set(handles.message_texts,'String','Updated (full-resolution) scan')
update_ptcld_plot(hObject,handles);
%--------------------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);

function update_ptcld_plot(hObject,handles)
MS     = str2double(get(handles.MStag,'String'));
%-------------------------------------------------------------------------------
% Perform translations of the points
if handles.cropped == 1
    [colobj,tdat,eclks] = update_PLY_pts(handles);
else
    colobj = handles.colobj;
    set(handles.message_texts,'String','Please Crop before rotating/translating');
end
%-------------------------------------------------------------------------------
% Update the figure
axes(handles.mainax);
azel = get(handles.mainax,'view');
cla
pcshow(colobj,'Markersize',MS);
hold on
if isfield(handles,'add_bestfit_els') == 1
    if handles.add_bestfit_els == 1
        set(handles.mainax,'view',azel);
        [Bout,elns,bestsclf] = try_bestfitting_of_nominal_elecs(handles);
        handles.Bout = Bout;
        handles.elns = elns;
        handles.bestsclf = bestsclf;
    end
end
%-------------------------------------------------------------------------------
axis equal
grid on
box on
FS = 12;
xlabel('X','fontsize',FS,'FontName','times')
ylabel('Y','fontsize',FS,'FontName','times')
zlabel('Z','fontsize',FS,'FontName','times')
set(gca,'FontSize',FS,'FontName','times')
set(handles.mainax,'view',azel);
% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
function ini_rx_Callback(hObject, eventdata, handles)
% hObject    handle to ini_rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ini_rx as text
%        str2double(get(hObject,'String')) returns contents of ini_rx as a double
update_ptcld_plot(hObject,handles)

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function ini_rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ini_rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
function ini_rz_Callback(hObject, eventdata, handles)
% hObject    handle to ini_rz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ini_rz as text
%        str2double(get(hObject,'String')) returns contents of ini_rz as a double
update_ptcld_plot(hObject,handles)

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function ini_rz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ini_rz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in set_ini_frm_btn.
function set_ini_frm_btn_Callback(hObject, eventdata, handles)
% hObject    handle to set_ini_frm_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%--------------------------------------------------------------------------
function zpl_x1_Callback(hObject, eventdata, handles)
% hObject    handle to zpl_x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpl_x1 as text
%        str2double(get(hObject,'String')) returns contents of zpl_x1 as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function zpl_x1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpl_x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
function zpl_x2_Callback(hObject, eventdata, handles)
% hObject    handle to zpl_x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpl_x2 as text
%        str2double(get(hObject,'String')) returns contents of zpl_x2 as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function zpl_x2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpl_x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
function zpl_y1_Callback(hObject, eventdata, handles)
% hObject    handle to zpl_y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpl_y1 as text
%        str2double(get(hObject,'String')) returns contents of zpl_y1 as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function zpl_y1_CreateFcn(hObject, eventdata, ~)
% hObject    handle to zpl_y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
function zpl_y2_Callback(hObject, eventdata, handles)
% hObject    handle to zpl_y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpl_y2 as text
%        str2double(get(hObject,'String')) returns contents of zpl_y2 as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function zpl_y2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpl_y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
function zpl_z_Callback(hObject, eventdata, handles)
% hObject    handle to zpl_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpl_z as text
%        str2double(get(hObject,'String')) returns contents of zpl_z as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function zpl_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpl_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 
% %--------------------------------------------------------------------------
% % --- Executes on button press in draw_pln_btn.
% function draw_pln_btn_Callback(hObject, eventdata, handles)
% % hObject    handle to draw_pln_btn (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% if handles.plotplane == 0
%     handles.plotplane = 1;
% elseif handles.plotplane == 1
%     handles.plotplane = 0;
% end
% update_ptcld_plot(hObject,handles)
% %-------------------------------------------------------------------------------
% % Update handles structure
% guidata(hObject, handles);


% --- Executes on button press in click_elecs.
% function click_elecs_Callback(hObject, eventdata, handles)
% % hObject    handle to click_elecs (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% %-------------------------------------------------------------------------------
% if handles.cropped == 1
%     
%     %-------------------------------------------------------------------------------
%     [colobj,tdat,eclks,cent0] = update_PLY_pts(handles);
%     vs   = colobj.Location;
%     cent = [ ...
%         str2double(get(handles.org_x,'String')) ...
%         str2double(get(handles.org_y,'String')) ...
%         str2double(get(handles.org_z,'String'))]';
%     %---------------------------------------------------------------------------
%     xbs = [str2double(get(handles.zpl_x1,'String')) ...
%         str2double(get(handles.zpl_x2,'String'))];
%     ybs = [str2double(get(handles.zpl_y1,'String')) ...
%         str2double(get(handles.zpl_y2,'String'))];
%     z0  = str2double(get(handles.zpl_z,'String'));
%     is  = find( vs(:,3) > z0);
%     %---------------------------------------------------------------------------
%     % set(gca,'handles.txt_msgs','Left Click to Add, Right Click to Remove, and Click out of axes to stop');
%     check = 0;
%     while check == 0
%         %-----------------------------------------------------------------------
%         % Click on a point
%         view(2)
%         [x,y,butt] = ginputWhite(1);
%         
%         if ( (x < xbs(1)) || (x > xbs(2)) || (y < ybs(1)) || (y > ybs(2)) )
%             check = 1;
%         elseif (butt == 1) || (butt == 3)
%             %-------------------------------------------------------------------
%             % Get the closest z-value
%             [tmp,i] = min(sqrt( (vs(is,1)-x).^2+(vs(is,2)-y).^2));
%             newp    = [x y vs(is(i),3)]';
%             %-------------------------------------------------------------------
%             % Convert the point back to the initial coordinates
%             el2  = str2double(get(handles.cur_el,'String'))*pi/180;
%             Rx   = [1 0 0; 0 cos(el2) -sin(el2); 0 sin(el2) cos(el2)];
%             newp = (Rx')*newp;
%             az2  = str2double(get(handles.cur_az,'String'))*pi/180;
%             Rz   = [cos(az2) sin(az2) 0;-sin(az2) cos(az2) 0; 0 0 1];
%             newp = (Rz')*newp;
%             newp = newp + cent;
%             %-----------------------------------------------------------------------
%             % Rotate the STL to the initial frame
%             az   = 0*pi/180; % str2double(get(handles.ini_rz,'String'))*pi/180;
%             el   = 90*pi/180; % str2double(get(handles.ini_rx,'String'))*pi/180;
%             Rz   = [cos(az) sin(az) 0;-sin(az) cos(az) 0; 0 0 1];
%             newp = (Rz')*newp;
%             Rx   = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
%             newp = (Rx')*newp;
%             %-------------------------------------------------------------------
%             newp = newp + cent0';
%             
%             if (butt == 1)
%                 %-------------------------------------------------------------------
%                 % Record the data and update the plot
%                 if isfield(handles,'eclks') == 1
%                     handles.eclks = [handles.eclks; newp'];
%                 else
%                     handles.eclks = newp';
%                 end
%             elseif (butt == 3)
%                 % Remove the closest electrode already clicked
%                 if size(handles.eclks,1) > 1
%                     [tmp,i0] = min( sum( (handles.eclks - repmat(newp',size(handles.eclks,1),1)).^2,2));
%                     handles.eclks(i0,:) = [];
%                 end
%             end
%             
%             %-------------------------------------------------------------------
%             % Update the plot
%             update_ptcld_plot(hObject,handles)
%             
%         else
%             error('stop: not coded')
%         end
%         
%     end
%     %---------------------------------------------------------------------------
%     % Update handles structure
%     guidata(hObject, handles);
% else
%     set(handles.message_texts,'String','You must Crop the PLY first!!')
% end


% --- Executes on slider movement.
function rx_slider_Callback(hObject, eventdata, handles)
% hObject    handle to rx_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
az = 360*get(handles.rx_slider,'Value');
el = 180*get(handles.rz_slider,'Value')-90;
% set(handles.mainax,'View',[az el])
set(handles.cur_az,'String',num2str(az))
update_ptcld_plot(hObject,handles)


% --- Executes during object creation, after setting all properties.
function rx_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function rz_slider_Callback(hObject, eventdata, handles)
% hObject    handle to rz_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

az = 360*get(handles.rx_slider,'Value');
el = 180*get(handles.rz_slider,'Value')-90;
% set(handles.mainax,'View',[az el])
set(handles.cur_el,'String',num2str(el))
update_ptcld_plot(hObject,handles)



% --- Executes during object creation, after setting all properties.
function rz_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rz_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function org_x_Callback(hObject, eventdata, handles)
% hObject    handle to org_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of org_x as text
%        str2double(get(hObject,'String')) returns contents of org_x as a double
update_ptcld_plot(hObject,handles)

% --- Executes during object creation, after setting all properties.
function org_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to org_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function org_y_Callback(hObject, eventdata, handles)
% hObject    handle to org_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of org_y as text
%        str2double(get(hObject,'String')) returns contents of org_y as a double
update_ptcld_plot(hObject,handles)

% --- Executes during object creation, after setting all properties.
function org_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to org_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function org_z_Callback(hObject, eventdata, handles)
% hObject    handle to org_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of org_z as text
%        str2double(get(hObject,'String')) returns contents of org_z as a double
update_ptcld_plot(hObject,handles)

% --- Executes during object creation, after setting all properties.
function org_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to org_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in save_btn.
function save_btn_Callback(hObject, eventdata, handles)
% hObject    handle to save_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%-------------------------------------------------------------------------------
%if it is a subscan, merge all the subscans to make a full point cloud
%before saving
if isfield(handles, 'is_sub')==1
if handles.is_sub==1;
    pc1_locs=handles.pc1.Location;
pc2_locs=handles.pc2.Location;
pc3_locs=handles.pc3.Location;
pc4_locs=handles.pc4.Location;

pc1_cols=handles.pc1.Color;
pc2_cols=handles.pc2.Color;
pc3_cols=handles.pc3.Color;
pc4_cols=handles.pc4.Color;

pc_full_locs=[pc1_locs; pc2_locs; pc3_locs; pc4_locs];
pc_full_cols=[pc1_cols; pc2_cols; pc3_cols; pc4_cols];
pc_full=pointCloud(pc_full_locs, 'Color', pc_full_cols);
handles.colobj=pc_full;
handles.is_sub=0;
%handles.scn_num=0;
  update_ptcld_plot(hObject,handles)
  set(handles.message_texts,'String', 'changed to fullview for saving')
guidata(hObject, handles);
end
end
%-------------------------------------------------------------------------------
colobj     = handles.colobj;
file_ini   = handles.file;
cropped    = handles.cropped;
%plotplane  = handles.plotplane;

if cropped == 0
    %---------------------------------------------------------------------------
    % If no cropping has been done, then nothing has really been done
    [file,path] = uiputfile([file_ini,'.mat'],'Save Cropped PLY Data');
    eval(['save ',path,'/',file,' file_ini colobj cropped plotplane'])
else
    %---------------------------------------------------------------------------
    % initial rotation frame
    ini_rx = 90*pi/180; % str2double(get(handles.ini_rx,'String'));
    ini_rz = 0; % str2double(get(handles.ini_rz,'String'));
    
    %---------------------------------------------------------------------------
    % If there is an origin defined subtract off the origin
    cent = [0 0 0];
    xbs  = [-0.1 0.1];
    ybs  = [-0.1 0.1];
    z0   = 0;
    %---------------------------------------------------------------------------
    tdat  = [];
    eclks = [];
    %targ_strs = get(handles.popupmenu1,'String');
    %---------------------------------------------------------------------------
    % If there is an origin defined subtract off the origin
    if isfield(handles,'tdat') == 1
        tdat      = handles.tdat;
        tdat(tdat==0)=NaN;
        if size(tdat, 1)<258
            tdat(end+1:258, :)=NaN; 
        end
    end
    if isfield(handles,'eclks') == 1
        eclks = handles.eclks;
    end
    if isfield(handles,'newtdats') == 1
        newtdats = handles.newtdats;
    else
        newtdats = [];
    end
    if isfield(handles,'azels') == 1
        azels = handles.azels;
    else
        azels = [];
    end
    if isfield(handles,'labtdats') == 1
        labtdats = handles.labtdats;
    else
        labtdats = [];
    end

    
    
    %---------------------------------------------------------------------------
    % If no cropping has been done, then nothing has really been done
    [file,path] = uiputfile([file_ini,'.mat'],'Save Cropped PLY Data');
    %------- HR Edit 10/15/21------------------------
       if isvarname(tdat)==0 % if the scan is cropped but no electrodes have been picked no electrodes have been clicked yet
  eval(['save ',path,file,' file_ini colobj cropped plotplane', ...
         ' eclks ini_rx ini_rz cent xbs ybs z0 azels'])
       end
     %-----------------------------------------------------
    eval(['save ',path,'/',file,' file_ini colobj cropped plotplane', ...
        ' tdat eclks ini_rx ini_rz cent xbs ybs z0 newtdats azels labtdats'])

end
%--------------------------------------------------------------------------
% --- Executes on button press in load_btn.
function load_btn_Callback(hObject, eventdata, handles)
% hObject    handle to load_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%addpath(genpath('S:/digihisto/Ethan/EKM_utility'))
addpath mfiles
set(handles.message_texts,'String','Loading Cropped PLY data')
%-------------------------------------------------------------------------------
[file,path] = uigetfile('*.mat','Load Segmentation Data');

dat = load(['',path,'/',file,'']);
set(handles.file_name_txtbx,'String',dat.file_ini)
handles.colobj    = dat.colobj;
handles.colobj0   = dat.colobj;
handles.file      = dat.file_ini;
handles.cropped   = dat.cropped;
handles.plotplane = dat.plotplane;
% Update handles structure
guidata(hObject, handles);

if handles.cropped == 1
    %----------------------------------------------------------------------
    % initial rotation frame
    % set(handles.ini_rx,'String',num2str(dat.ini_rx));
    % set(handles.ini_rz,'String',num2str(dat.ini_rz));
    handles.ini_rx = 90*pi/180;
    handles.ini_rz = 0;
    
    %----------------------------------------------------------------------
    % If there is an origin defined subtract off the origin
    % set(handles.org_x,'String',num2str(dat.cent(1)));
    % set(handles.org_y,'String',num2str(dat.cent(2)));
    % set(handles.org_z,'String',num2str(dat.cent(3)));
    %     if handles.plotplane == 1
    %         set(handles.zpl_x1,'String',num2str(dat.xbs(1)));
    %         set(handles.zpl_x2,'String',num2str(dat.xbs(2)));
    %         set(handles.zpl_y1,'String',num2str(dat.ybs(1)));
    %         set(handles.zpl_y2,'String',num2str(dat.ybs(2)));
    %         set(handles.zpl_z,'String',num2str(dat.z0));
    %     end
    
    %---------------------------------------------------------------------------
    if isfield(handles,'eclks')
    if size(dat.eclks,1) > 0
        handles.eclks = dat.eclks;
    else
        
            handles = rmfield(handles,'eclks');
    end
    end
    % dat.tdat = [];
    if isfield(dat, 'tdat')==1
    if size(dat.tdat,1) > 0
        %targ_strs = get(handles.popupmenu1,'String');
        tdat      = NaN*ones(261,3);
        tdat(1:size(dat.tdat,1),:) = dat.tdat;
        handles.tdat = tdat;
        pc_in=handles.colobj;
        [pc_out]=make_targets_green(pc_in, tdat, 25);
        handles.colobj=pc_out;
    else
        if isfield(handles,'tdat')
            handles = rmfield(handles,'tdat');
        end
    end
    end
    if isfield(dat,'newtdats') == 1
        if size(dat.newtdats,1) > 0
            handles.newtdats = dat.newtdats;
            newtdats=handles.newtdats;
            pc_in=handles.colobj;
           [pc_out]=make_elecs_red(pc_in, newtdats, 25);
           handles.colobj=pc_out;
        end
    end
    if isfield(dat,'azels') == 1
        handles.azels = dat.azels;        
    end
    
    %     if isfield(dat,'targ_strs') == 1
    %         set(handles.popupmenu1,'String',dat.targ_strs);
    %     end
end
%-------------------------------------------------------------------------------
% Update the figure
update_ptcld_plot(hObject,handles)

%-------------------------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);



% % --------------------------------------------------------------------
% function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% % hObject    handle to uitoggletool1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% azel = get(handles.mainax,'View');
% set(handles.cur_az,'String',num2str(azel(1)))
% set(handles.cur_el,'String',num2str(azel(2)))


%--------------------------------------------------------------------------
function cur_az_Callback(hObject, eventdata, handles)
% hObject    handle to cur_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cur_az as text
%        str2double(get(hObject,'String')) returns contents of cur_az as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function cur_az_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cur_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
function cur_el_Callback(hObject, eventdata, handles)
% hObject    handle to cur_el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cur_el as text
%        str2double(get(hObject,'String')) returns contents of cur_el as a double

% %--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function cur_el_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cur_el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in view2Dtag.
% function view2Dtag_Callback(hObject, eventdata, handles)
% % hObject    handle to view2Dtag (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% view(2)
% %--------------------------------------------------------------------------
% % --- Executes on button press in view3Dtag.
% function view3Dtag_Callback(hObject, eventdata, handles)
% % hObject    handle to view3Dtag (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% view(3)
%--------------------------------------------------------------------------


function MStag_Callback(hObject, eventdata, handles)
% hObject    handle to MStag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MStag as text
%        str2double(get(hObject,'String')) returns contents of MStag as a double
update_ptcld_plot(hObject,handles)

% --- Executes during object creation, after setting all properties.
function MStag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MStag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in set_targloc.
% function set_targloc_Callback(hObject, eventdata, handles)
% % hObject    handle to set_targloc (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% %-------------------------------------------------------------------------------
% if handles.cropped == 1
%     
%     %----------------------------------------------------------------------
%     [colobj,tdat_tmp,eclks_tmp,cent0] = update_PLY_pts(handles);
%     vs   = colobj.Location;
%     cent = [0 0 0]';
%     xbs  = [-0.1 0.1];
%     ybs  = [-0.1 0.1];
%     z0   = 0;
%     is  = find( vs(:,3) > z0);
%     %---------------------------------------------------------------------------
%     targ_strs = get(handles.popupmenu1,'String');
%     i_t       = get(handles.popupmenu1,'Value');
%     set(handles.message_texts,'String',['Click on ',targ_strs{i_t}]);
%     %-----------------------------------------------------------------------
%     % Click on a point
%     view(2)
%     [x,y,butt] = ginputWhite(1);
%     
%     if ( (x < xbs(1)) || (x > xbs(2)) || (y < ybs(1)) || (y > ybs(2)) )
%         % Do nothing if they click outside of the axis
%     elseif (butt == 3)
%         if isfield(handles,'tdat') == 1
%             handles.tdat(i_t,:) = NaN;
%         else
%             handles.tdat = NaN*ones(length(targ_strs),3);
%         end
%     elseif (butt == 1)
%         %-------------------------------------------------------------------
%         % Get the closest z-value
%         [tmp,i] = min(sqrt( (vs(is,1)-x).^2+(vs(is,2)-y).^2));
%         newp    = [x y vs(is(i),3)]';
%         %-------------------------------------------------------------------
%         % Convert the point back to the initial coordinates
%         el2  = str2double(get(handles.cur_el,'String'))*pi/180;
%         Rx   = [1 0 0; 0 cos(el2) -sin(el2); 0 sin(el2) cos(el2)];
%         newp = (Rx')*newp;
%         az2  = str2double(get(handles.cur_az,'String'))*pi/180;
%         Rz   = [cos(az2) sin(az2) 0;-sin(az2) cos(az2) 0; 0 0 1];
%         newp = (Rz')*newp;
%         newp = newp + cent;
%         %-----------------------------------------------------------------------
%         % Rotate the STL to the initial frame
%         az   = handles.ini_rz;
%         el   = handles.ini_rx;
%         Rz   = [cos(az) sin(az) 0;-sin(az) cos(az) 0; 0 0 1];
%         newp = (Rz')*newp;
%         Rx   = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
%         newp = (Rx')*newp;
%         %-------------------------------------------------------------------
%         newp = newp + cent0';
%         %-------------------------------------------------------------------
%         if isfield(handles,'tdat') == 1
%             handles.tdat(i_t,:) = newp;
%         else
%             handles.tdat = NaN*ones(length(targ_strs),3);
%             handles.tdat(i_t,:) = newp;
%         end
%         
%         %-------------------------------------------------------------------
%         % Update the plot
%         update_ptcld_plot(hObject,handles)
%         
%     else
%         error('stop: not coded')
%     end
%     
%     
%     %---------------------------------------------------------------------------
%     % Update handles structure
%     guidata(hObject, handles);
% else
%     set(handles.message_texts,'String','You must Crop the PLY first!!')
% end


% % --- Executes on button press in load_new_targlist.
% function load_new_targlist_Callback(hObject, eventdata, handles)
% % hObject    handle to load_new_targlist (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% %-------------------------------------------------------------------------------
% % Load a new list file
% [file,path] = uigetfile('*.txt');
% 
% if isequal(file,0) || isequal(path,0)
%     set(handles.message_texts,'String','User pressed cancel')
% else
%     %---------------------------------------------------------------------------
%     dat = strsplit(fileread([path,'/',file]));
%     set(handles.popupmenu1,'String',dat)
%     set(handles.message_texts,'String','List of Possible Targets Updated')
% end


% % --- Executes on button press in rmv_seltarg.
% function rmv_seltarg_Callback(hObject, eventdata, handles)
% % hObject    handle to rmv_seltarg (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% %-------------------------------------------------------------------------------
% i = get(handles.popupmenu1,'Value');
% handles.tdat(i,:) = NaN;
% guidata(hObject, handles);


%--- Executes on button press in addbest_fit_els.
function addbest_fit_els_Callback(hObject, eventdata, handles)
% hObject    handle to addbest_fit_els (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % if isfield(handles,'add_bestfit_els') == 0
handles.add_bestfit_els = 0;
% load dat/UCL_head_mesh_256eeg_elecs msh el256 eeg_mps
nomdat = load_nominal_cap(0);
handles.msh     = nomdat.msh;
handles.el256   = nomdat.el256;
handles.eeg_mps = nomdat.eeg_mps(1:257,:); % convert from cm to mm

%     figure
%     plot3(handles.eeg_mps(:,1),handles.eeg_mps(:,2),handles.eeg_mps(:,3),'.k','markersize',12)
% % end

if handles.add_bestfit_els == 1
    handles.add_bestfit_els = 0;
elseif handles.add_bestfit_els == 0
    handles.add_bestfit_els = 1;
end
%-------------------------------------------------------------------------------
% Update the plot
update_ptcld_plot(hObject,handles)

%-------------------------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);


% % --- Executes on button press in run_prj_img_proc.
% function run_prj_img_proc_Callback(hObject, eventdata, handles)
% % hObject    handle to run_prj_img_proc (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% 
% if isfield(handles,'add_bestfit_els') == 1
%     if handles.add_bestfit_els == 1
%         %-----------------------------------------------------------------------
%         % Perform translations of the points
%         [colobj,tdat,~,cent0] = update_PLY_pts(handles);
%         
%         %-----------------------------------------------------------------------
%         % Get the indices of the located target points
%         is_fnd = find( (isnan(tdat(:,1)) == 0) );
%         is_fnd = is_fnd( is_fnd < 257);
%         
%         if length(is_fnd) >= 3
%             %-------------------------------------------------------------------
%             % Run the projection algorithm
%             %             save testdat handles
%             %             clear
%             % % newtdats = run_prj_alg(colobj,[], 0);
%             
%             %incorrect input args here? 
%             [newtdats] = run_multiview_prj_alg_v5(handles);
%             
%             handles.newtdats = newtdats;
%             handles.eclks    = newtdats;
%             handles.azels    = [];
%             update_ptcld_plot(hObject,handles)
%             
%             %-------------------------------------------------------------------
%             % Update handles structure
%             guidata(hObject, handles);
%             %             set(handles.message_texts,'String','For each new potential target point click on the corresponding electrode number or right click if its wrong')
%             %             vs       = colobj.Location;
%             %             cent = [ ...
%             %                 str2double(get(handles.org_x,'String')) ...
%             %                 str2double(get(handles.org_y,'String')) ...
%             %                 str2double(get(handles.org_z,'String'))]';
%             %             %---------------------------------------------------------------------------
%             %             xbs = [str2double(get(handles.zpl_x1,'String')) ...
%             %                 str2double(get(handles.zpl_x2,'String'))];
%             %             ybs = [str2double(get(handles.zpl_y1,'String')) ...
%             %                 str2double(get(handles.zpl_y2,'String'))];
%             %             z0  = str2double(get(handles.zpl_z,'String'));
%             %             is  = find( vs(:,3) > z0);
%             %
%             %             for n = 1:size(newtdats,1)
%             %                 % Skip the work if the new possible target point is within
%             %                 % 5 mm of a current target point
%             %                 minds = min(sqrt(sum( (tdat - newtdats(n,:)).^2)),[],'omitnan');
%             %                 if minds*1000 < 5
%             %                     set(handles.message_texts,'String',['New target ',num2str(n),' of ',num2str(size(newtdats,1)),': Already found'])
%             %                 else
%             %                     %-----------------------------------------------------------
%             %                     % Plot the current new target point
%             %                     axes(handles.mainax);
%             %                     hold on
%             %                     fac = 1.05;
%             %                     plot3(fac*newtdats(n,1),fac*newtdats(n,2),fac*newtdats(n,3),'.g','markersize',28)
%             %                     istmp = setdiff(1:size(newtdats,1),n);
%             %                     plot3(fac*newtdats(istmp,1),fac*newtdats(istmp,2),fac*newtdats(istmp,3),'.m','markersize',12)
%             %
%             %
%             %
%             %                     %-----------------------------------------------------------
%             %                     % Click on a point
%             %                     view(2)
%             %                     [x,y,butt] = ginputWhite(1);
%             %
%             %                     if ( (x < xbs(1)) || (x > xbs(2)) || (y < ybs(1)) || (y > ybs(2)) ) || (butt == 3)
%             %                         check = 1;
%             %                     elseif (butt == 1)
%             %                         %-------------------------------------------------------
%             %                         % Get the closest z-value
%             %                         [tmp,i] = min(sqrt( (vs(is,1)-x).^2+(vs(is,2)-y).^2));
%             %                         newp0   = [x y vs(is(i),3)]';
%             %                         newp    = newtdats(n,:)';
%             %                         %-------------------------------------------------------
%             %                         % Convert the point back to the initial coordinates
%             %                         el2  = str2double(get(handles.cur_el,'String'))*pi/180;
%             %                         Rx   = [1 0 0; 0 cos(el2) -sin(el2); 0 sin(el2) cos(el2)];
%             %                         newp = (Rx')*newp;
%             %                         az2  = str2double(get(handles.cur_az,'String'))*pi/180;
%             %                         Rz   = [cos(az2) sin(az2) 0;-sin(az2) cos(az2) 0; 0 0 1];
%             %                         newp = (Rz')*newp;
%             %                         newp = newp + cent;
%             %                         %-------------------------------------------------------
%             %                         % Rotate the STL to the initial frame
%             %                         az   = str2double(get(handles.ini_rz,'String'))*pi/180;
%             %                         el   = str2double(get(handles.ini_rx,'String'))*pi/180;
%             %                         Rz   = [cos(az) sin(az) 0;-sin(az) cos(az) 0; 0 0 1];
%             %                         newp = (Rz')*newp;
%             %                         Rx   = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
%             %                         newp = (Rx')*newp;
%             %                         %-------------------------------------------------------
%             %                         newp = newp + cent0';
%             %                         %-------------------------------------------------------
%             %                         % Find the closest electrode number that you picked on
%             %                         % It should be the closest point to the nominal,
%             %                         % best fit electrodes
%             %                         [Bout,elns] = try_bestfitting_of_nominal_elecs(handles);
%             %
%             %                         %                         Bout     = handles.Bout;
%             %                         %                         elns     = handles.elns;
%             %                         istmp2   = find(Bout(:,3) > 0);
%             %                         Bout     = Bout(istmp2,:);
%             %                         elns     = elns(istmp2);
%             %
%             %
%             %                         [tmp,i0] = min(sqrt(sum( (Bout(:,1:3) - repmat(double(newp0'),size(Bout,1),1)).^2,2)));
%             %                         Bout(i0,1:3)*1000
%             %                         newp0'*1000
%             %                         tmp*1000
%             %                         elns(i0)
%             %                         handles.tdat(elns(i0),:) = newp;
%             %
%             %
%             %                         %---------------------------------------------------------------
%             %                         % Update the plot
%             %                         update_ptcld_plot(hObject,handles)
%             %                     end
%             %                     newtdats(n,:) = NaN;
%             %
%             %                 end
%             % end
%         else
%             set(handles.message_texts,'String','Pick the manual labeled electrodes')
%             
%         end    
% 
%     else
%         set(handles.message_texts,'String','Make sure to add best electrode fits')
%     end
% else
%     set(handles.message_texts,'String','Make sure to add best electrode fits')        
% end


% --- Executes on button press in single_run_prj_img_proc.
function single_run_prj_img_proc_Callback(hObject, eventdata, handles)
% hObject    handle to single_run_prj_img_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



        % Run the projection algorithm
        %             save testdat
        %             clear
        % newtdats = run_prj_alg(colobj);
        if isfield(handles,'newtdats') == 1        
            newtdats = handles.newtdats;
        else
            newtdats = [];
        end
        [newtdats,colobj_out] = run_single_prj_alg_v2_sub(newtdats,handles, handles.colobj);
        handles.colobj=colobj_out;
        handles.newtdats = newtdats;
        handles.eclks    = newtdats;
        update_ptcld_plot(hObject,handles)
        
            
        %if you are working with a subscan, update the corresponding subscan field
        if isfield(handles, 'is_sub')==1
        if handles.is_sub==1
            scn_num=handles.scn_num;
            switch scn_num
                case 1
                    handles.pc1=handles.colobj;
                   [pc2]=update_pc_intersection(hObject, handles, 1, 2);
                   [pc3]=update_pc_intersection(hObject, handles, 1, 3) ;
                   [pc4]=update_pc_intersection(hObject, handles, 1, 4) ;
                   handles.pc2=pc2;
                   handles.pc3=pc3;
                   handles.pc4=pc4;
              case 2
                    handles.pc2=handles.colobj;
                   [pc1]=update_pc_intersection(hObject, handles, 2, 1);
                   [pc3]=update_pc_intersection(hObject, handles, 2, 3) ;
                   [pc4]=update_pc_intersection(hObject, handles, 2, 4) ;
                   handles.pc1=pc1;
                   handles.pc3=pc3;
                   handles.pc4=pc4;
                    
              case 3
                    handles.pc3=handles.colobj;
                   [pc1]=update_pc_intersection(hObject, handles, 3, 1);
                   [pc2]=update_pc_intersection(hObject, handles, 3, 2) ;
                   [pc4]=update_pc_intersection(hObject, handles, 3, 4) ;
                   handles.pc1=pc1;
                   handles.pc2=pc2;
                   handles.pc4=pc4;
              case 4
                    handles.pc4=handles.colobj;
                   [pc1]=update_pc_intersection(hObject, handles, 4, 1);
                   [pc2]=update_pc_intersection(hObject, handles, 4, 2) ;
                   [pc3]=update_pc_intersection(hObject, handles, 4, 3) ;
                   handles.pc1=pc1;
                   handles.pc2=pc2;
                   handles.pc3=pc3;
            end
        end
      
        
        end
        %-------------------------------------------------------------------
        % Update handles structure
        guidata(hObject, handles);



% % --- Executes on button press in multilabels.
% function multilabels_Callback(hObject, eventdata, handles)
% % hObject    handle to multilabels (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% check = 0;
% if isfield(handles,'add_bestfit_els') == 1
%     if handles.add_bestfit_els == 1
%         if isfield(handles,'newtdats') == 1
%             if size(handles.newtdats,1) > 0
%                 check = 1;
%             end
%         end
%     end
% end
% 
% if check == 1
%     %---------------------------------------------------------------------------
%     % Run the projection labeling algorithm
%     
%     [labtdats,tdats] = run_multiview_labeling_prj_alg(handles);
%     handles.labtdats = labtdats;
%     handles.tdat     = tdats;
%         
%     update_ptcld_plot(hObject,handles)
%     
%     %---------------------------------------------------------------------------
%     % Update handles structure
%     guidata(hObject, handles);
% else
%     set(handles.message_texts,'String','Make sure electrodes were previously labeled and a best-fit of nominal electrodes was done')
%     
% end

function figure1_KeyPressFcn(hObject, eventdata, handles)
disp('here')
uiresume;


% --- Executes on button press in singledir_remove.
function singledir_remove_Callback(hObject, eventdata, handles)
% hObject    handle to singledir_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to single_run_prj_img_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if (isfield(handles,'newtdats') == 1)
        %-------------------------------------------------------------------
        % Run the projection algorithm
        %             save testdat
        %             clear
        % newtdats = run_prj_alg(colobj);
        newtdats = handles.newtdats;
        [newtdats_out, colobj_out] = run_single_prj_removals(newtdats,handles);
        
        handles.newtdats = newtdats_out;
        handles.eclks    = newtdats;
         handles.colobj=colobj_out;
        update_ptcld_plot(hObject,handles)
       
        %-------------------------------------------------------------------
        % Update handles structure
        guidata(hObject, handles);
    end



% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function figure1_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in find_missing_els.
function find_missing_els_Callback(hObject, eventdata, handles)
% hObject    handle to find_missing_els (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

check = 0;
if isfield(handles,'add_bestfit_els') == 1
    if handles.add_bestfit_els == 1
        
        check = 1;
        
    end
end

if check == 1
    %---------------------------------------------------------------------------
    % Run the projection labeling algorithm
    
    [labtdats,tdats] = run_missing_els_prj_alg(handles);
    handles.labtdats = labtdats;
    handles.tdat     = tdats;

    update_ptcld_plot(hObject,handles)
    
    %---------------------------------------------------------------------------
    % Update handles structure
    guidata(hObject, handles);
else
    set(handles.message_texts,'String','Make sure electrodes were previously labeled and a best-fit of nominal electrodes was done')
    
end

%--------------------------------------------------------------------------
% --- Executes on button press in deidentify_face.
function deidentify_face_Callback(hObject, eventdata, handles)
% hObject    handle to deidentify_face (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------
% De-identify the face
set(handles.message_texts,'String','Make sure the view is centered over the identifyable features (face, eyes, nose, etc.)')
%--------------------------------------------------------------------------
handles = run_deidentify_alg(handles);
update_ptcld_plot(hObject,handles)

%-------------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in singlelabels.
% function singlelabels_Callback(hObject, eventdata, handles)
% % hObject    handle to singlelabels (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% check = 0;
% if isfield(handles,'add_bestfit_els') == 1
%     if handles.add_bestfit_els == 1
%         if isfield(handles,'newtdats') == 1
%             if size(handles.newtdats,1) > 0
%                 check = 1;
%             end
%         end
%     end
% end
% 
% if check == 1
%     %---------------------------------------------------------------------------
%     % Run the projection labeling algorithm
%     
%     [labtdats,tdats] = run_singleiview_labeling_prj_alg(handles);
%     handles.labtdats = labtdats;
%     handles.tdat     = tdats;
% 
%     update_ptcld_plot(hObject,handles)
%     
%     %---------------------------------------------------------------------------
%     % Update handles structure
%     guidata(hObject, handles);
% else
%     set(handles.message_texts,'String','Make sure electrodes were previously labeled and a best-fit of nominal electrodes was done')
%     
% end

% --- Executes on button press in autolabeling.
% function autolabeling_Callback(hObject, eventdata, handles)
% % hObject    handle to autolabeling (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% check = 0;
% if isfield(handles,'add_bestfit_els') == 1
%     if handles.add_bestfit_els == 1
%         if isfield(handles,'newtdats') == 1
%             if size(handles.newtdats,1) > 0
%                 check = 1;
%             end
%         end
%     end
% end
% 
% if check == 1
%     %----------------------------------------------------------------------
%     % Run the projection labeling algorithm
%     [tdat,elunlab] = run_autlabel_funcs_for_GUI(handles);
%     %----
%     % handles.tdat      = tdat;        
%     % handles.newtdats  = elunlab;
%     % handles.eclks     = elunlab;
%     %----
%     elfound           = find(isnan(tdat(:,1))==0);
%     newtdats          = [tdat(elfound,:); elunlab];
%     handles.newtdats  = newtdats;
%     handles.eclks     = newtdats;
%     handles.labtdats  = [elfound; NaN*elunlab(:,1)];
%     % labtdats          = NaN*tdat(:,1);
%     % elfound           = find(isnan(tdat(:,1))==0);
%     % labtdats(elfound) = elfound;
%     % handles.labtdats  = labtdats;        
%     % handles.labtdats  = find(isnan(tdat(:,1))==0);   
%     update_ptcld_plot(hObject,handles)
%     
%     %----------------------------------------------------------------------
%     % Update handles structure
%     guidata(hObject, handles);
% else
%     set(handles.message_texts,'String','Make sure electrodes were previously labeled and a best-fit of nominal electrodes was done')
%     
% end


% --- Executes on button press in subdivide_btn.
function subdivide_btn_Callback(hObject, eventdata, handles)
% hObject    handle to subdivide_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n_part=4;
handles.colobj_orig=handles.colobj;
[colobj,tdat,eclks,cent0] = update_PLY_pts(handles);

[pc_sub, proj_data]= project_scan_onto_axes_and_subdivide_v4(colobj, n_part) ;
locs=handles.colobj0.Location;
for k=1:n_part
    xyz=locs(pc_sub{k}, :);
    pc_test{k}=pointCloud(xyz, 'Color', colobj.Color(pc_sub{k}, :));
end
disp('got here')


handles.pc_sub1=pc_sub{1};
handles.pc_sub2=pc_sub{2};
handles.pc_sub3=pc_sub{3};
handles.pc_sub4=pc_sub{4};

handles.pc1=pc_test{1};
handles.pc2=pc_test{2};
handles.pc3=pc_test{3};
handles.pc4=pc_test{4};
%handles.proj_data=proj_data;
set(handles.message_texts,'String', 'scans are subdivided, choose a scan to process')
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in scn1.
function scn1_Callback(hObject, eventdata, handles)
% hObject    handle to scn1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.colobj=handles.pc1;
handles.is_sub=1;
handles.scn_num=1;
 update_ptcld_plot(hObject,handles)
 set(handles.message_texts,'String', 'scan 1')
 update_small_axis(hObject, handles);
guidata(hObject, handles);


% --- Executes on button press in scn2.
function scn2_Callback(hObject, eventdata, handles)
% hObject    handle to scn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.colobj=handles.pc2;
handles.is_sub=1;
handles.scn_num=2;
update_ptcld_plot(hObject,handles)
set(handles.message_texts,'String', 'scan 2')
update_small_axis(hObject, handles);
guidata(hObject, handles);


% --- Executes on button press in scn3.
function scn3_Callback(hObject, eventdata, handles)
% hObject    handle to scn3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.colobj=handles.pc3;
handles.is_sub=1;
handles.scn_num=3;
update_ptcld_plot(hObject,handles)
set(handles.message_texts,'String', 'scan 3')
update_small_axis(hObject, handles);
guidata(hObject, handles);


% --- Executes on button press in scn4.
function scn4_Callback(hObject, eventdata, handles)
% hObject    handle to scn4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.colobj=handles.pc4;
handles.is_sub=1;
handles.scn_num=4;
update_ptcld_plot(hObject,handles)
set(handles.message_texts,'String', 'scan 4')
update_small_axis(hObject, handles);
guidata(hObject, handles);


% --- Executes on button press in fullview.
function fullview_Callback(hObject, eventdata, handles)
% hObject    handle to fullview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pc1_locs=handles.pc1.Location;
pc2_locs=handles.pc2.Location;
pc3_locs=handles.pc3.Location;
pc4_locs=handles.pc4.Location;

pc1_cols=handles.pc1.Color;
pc2_cols=handles.pc2.Color;
pc3_cols=handles.pc3.Color;
pc4_cols=handles.pc4.Color;

pc_full_locs=[pc1_locs; pc2_locs; pc3_locs; pc4_locs];
pc_full_cols=[pc1_cols; pc2_cols; pc3_cols; pc4_cols];
pc_full=pointCloud(pc_full_locs, 'Color', pc_full_cols);
handles.colobj=pc_full;
handles.is_sub=0;

update_ptcld_plot(hObject,handles)
update_small_axis(hObject, handles)
set(handles.message_texts,'String', 'Full View')
guidata(hObject, handles);


% --- Executes on button press in single_targ_label.
function single_targ_label_Callback(hObject, eventdata, handles)
% hObject    handle to single_targ_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
         axes(handles.mainax);
         im=flip(imread('EGI_black_electrode_map.jpg'));
         imagesc(im)
         set(gca, 'YDir','normal')
          set(gca, 'XDir','normal')
          set(handles.message_texts, 'Visible', 'off')
         %title('reference map (black electrodes in red)')
      if isfield(handles, 'is_sub')==1
          if handles.is_sub==1
         update_small_axis(hObject, handles);
          end
      end
          %pcshow(colobj,'Markersize',MS);
        % Run the projection algorithm
% %         %             save testdat
        %             clear
        % newtdats = run_prj_alg(colobj);
        if isfield(handles,'tdat') == 1        
            tdat = handles.tdat;
        else
            tdat = [];
        end
        [tdat_out, colobj_out] = run_single_targ_alg_v2_sub(tdat,handles, handles.colobj);
        set(handles.message_texts, 'Visible', 'on')
        
        handles.tdat = tdat_out;
        handles.colobj=colobj_out;
        %handles.eclks    = newtdats;
        update_ptcld_plot(hObject,handles)
        
        %if you are working with a subscan, update the corresponding subscan field
        if isfield(handles, 'is_sub')==1
        if handles.is_sub==1
            scn_num=handles.scn_num;
            switch scn_num
                case 1
                    handles.pc1=handles.colobj;
                   [pc2]=update_pc_intersection(hObject, handles, 1, 2);
                   [pc3]=update_pc_intersection(hObject, handles, 1, 3) ;
                   [pc4]=update_pc_intersection(hObject, handles, 1, 4) ;
                   handles.pc2=pc2;
                   handles.pc3=pc3;
                   handles.pc4=pc4;
              case 2
                    handles.pc2=handles.colobj;
                   [pc1]=update_pc_intersection(hObject, handles, 2, 1);
                   [pc3]=update_pc_intersection(hObject, handles, 2, 3) ;
                   [pc4]=update_pc_intersection(hObject, handles, 2, 4) ;
                   handles.pc1=pc1;
                   handles.pc3=pc3;
                   handles.pc4=pc4;
                    
              case 3
                    handles.pc3=handles.colobj;
                   [pc1]=update_pc_intersection(hObject, handles, 3, 1);
                   [pc2]=update_pc_intersection(hObject, handles, 3, 2) ;
                   [pc4]=update_pc_intersection(hObject, handles, 3, 4) ;
                   handles.pc1=pc1;
                   handles.pc2=pc2;
                   handles.pc4=pc4;
              case 4
                    handles.pc4=handles.colobj;
                   [pc1]=update_pc_intersection(hObject, handles, 4, 1);
                   [pc2]=update_pc_intersection(hObject, handles, 4, 2) ;
                   [pc3]=update_pc_intersection(hObject, handles, 4, 3) ;
                   handles.pc1=pc1;
                   handles.pc2=pc2;
                   handles.pc3=pc3;
            end
        end
      
        
        end
        %-------------------------------------------------------------------
        % Update handles structure
        guidata(hObject, handles);


% --- Executes on button press in deidentify_face.
function pushbutton33_Callback(hObject, ~, ~)
% hObject    handle to deidentify_face (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in merge_btn.
function merge_btn_Callback(hObject, eventdata, handles)
% hObject    handle to merge_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get path and filename info:
if isfield(handles, 'subj_name')==0
set(handles.message_texts, 'String', 'type subject name first')
    set(handles.message_texts, 'String', 'Type subject name')
    set(handles.subj_name_btn, 'Visible', 'on')
    uiwait;
    set(handles.subj_name_btn, 'Visible', 'off')
    subj_name=get(handles.subj_name_btn,'String');
else
    subj_name=handles.subj_name_btn;
end
base_pth=[cd,'\dat\', subj_name]
scnpth_raw =[base_pth, '\', subj_name,'_raw'];
cd(scnpth_raw)
file_info=dir;
n_files=size(file_info, 1);
n_files=n_files-2;
for i=1:n_files
    scnfls_raw{i}=file_info(i+2).name;
end

scnpth  =[base_pth,'\', subj_name,'_elfound'];
for i=1:max(size(scnfls_raw))
    temp_name=scnfls_raw{i};
    temp_name=temp_name(1:end-4);
    temp_name=[temp_name,'.mat'];
    scnfls{i}=temp_name;
end   
%path to cropped .mat files
scnpth_cropped=[base_pth,'\', subj_name, '\_cropped'];
%set(handles.message_texts,'String','Extracting Triangles... please wait')
%extract_tris_for_gui(scnpth, scnpth_raw, scnfls, scnfls_raw)
%get electrode thickness
set(handles.elthick_txt, 'Visible', 'on')
set(handles.message_texts,'String','Type electrode thickness (mm)');
uiwait
elthick=str2num(get(handles.elthick_txt, 'string'));
uiresume;
set(handles.message_texts,'String','Merging Scans... please wait')
%update handles
guidata(hObject, handles)
%see what is in the reference folder
valpth = [base_pth,'\', subj_name, '_pgdat'];
%make a saving folder if it doesnt already exist
if isfolder([base_pth '\', subj_name, '_saved'])==0
mkdir([base_pth, '\', subj_name, '_saved'])
end
savepth=[base_pth, '\', subj_name, '_saved'];
if isfolder(valpth)==0 %if no reference data is available
%proceed with merging
[colobj_mrg, tdat_mrg, allscans]= merge_after_labeling_for_gui_vf_no_pg(scnpth, scnfls, scnpth_raw, scnfls_raw, subj_name, elthick)
%see if photogrammetry data/comparison is needed  
handles.colobj=colobj_mrg;
handles.tdat=tdat_mrg;
set(handles.message_texts, 'String', 'done merging')
guidata(hObject, handles)
else
  tmp=dir(valpth);
  name=tmp(3).name;
  file_ext=name(end-3:end)
%merge
[colobj_mrg, tdat_mrg, allscans]= merge_after_labeling_for_gui_vf_no_pg(scnpth, scnfls, scnpth_raw, scnfls_raw, subj_name, elthick, savepth)
 set(handles.message_texts, 'String', 'done merging')
 guidata(hObject, handles)
 
%compare with photogrammetry. Handles file types .mat and .sfp 
[rms_err, rms_full]= comp_pg_tdat_vf(tdat_mrg, colobj_mrg, valpth, subj_name, file_ext, savepth)
    set(handles.message_texts, 'String', 'done comparing reference')
     guidata(hObject, handles)
  end

    guidata(hObject, handles)
%clear all figures that arent the GUI 
isfig=1;
c=2;
while isfig==1
    if isgraphics(c)==1
close(c)
c=c+1;
    else
        isfig=0;
    end
end

function subj_name_btn_Callback(hObject, eventdata, handles)
% hObject    handle to subj_name_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subj_name_btn as text
%        str2double(get(hObject,'String')) returns contents of subj_name_btn as a double
  subj_name=get(hObject,'String');
  set(handles.message_texts,'String',['subject name is: ' subj_name]);
  handles.subj_name=subj_name;
   guidata(hObject, handles)
  uiresume;
 
%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function subj_name_btn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subj_name_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on selection change in downsample_list.
function downsample_list_Callback(hObject, eventdata, handles)
% hObject    handle to downsample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dat=cellstr(get(hObject,'String'))
select=get(hObject, 'Value');
down_fac=str2num(cell2mat(dat(select)))
colobj=handles.colobj;
 [colobj,~,eclks] = update_PLY_pts(handles);
colobj_ds=pcdownsample(colobj, 'random', down_fac);
handles.colobj=colobj_ds;
handles.colobj0=colobj_ds;

pcshow(colobj)
view(2)

set(handles.message_texts, 'String', 'downsampled point cloud')
update_ptcld_plot(hObject,handles)
guidata(hObject, handles);


%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function downsample_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to downsample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
% --- Executes on button press in clear_tdat_btn.
function clear_tdat_btn_Callback(hObject, eventdata, handles)
% hObject    handle to clear_tdat_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tdat=handles.tdat;
pc_in=handles.colobj;
[pc_out]=make_old_targets_black(pc_in, tdat, 25);
handles.colobj=pc_out;
% newtdats=handles.newtdats;
% %eclks=handles.eclks;
% eclks=[];
% newtdats=[];
% handles.newtdats=newtdats;
% handles.eclks=eclks;
tdat(1:end, :)=NaN;
handles.tdat=tdat;
set(handles.message_texts, 'String', 'tdats cleared')
update_ptcld_plot(hObject,handles)
% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% --- Executes on button press in double_check.
function double_check_Callback(hObject, eventdata, handles)
% hObject    handle to double_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
hold on
colobj=handles.colobj;
tdat=handles.tdat;
newtdats=handles.newtdats;
pcshow(colobj)
plot3(tdat(:, 1), tdat(:, 2), tdat(:, 3), 'r.', 'MarkerSize', 30)
plot3(newtdats(:, 1), newtdats(:, 2), newtdats(:, 3), 'b.', 'MarkerSize', 20)

%--------------------------------------------------------------------------
% --- Executes on button press in clear_els.
function clear_els_Callback(hObject, eventdata, handles)
% hObject    handle to clear_els (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to clear_tdat_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newtdats=handles.newtdats;
pc_in=handles.colobj;
[pc_out]=make_old_els_black(pc_in);
handles.colobj=pc_out;
% newtdats=handles.newtdats;
% %eclks=handles.eclks;
% eclks=[];
% newtdats=[];
% handles.newtdats=newtdats;
% handles.eclks=eclks;
newtdats(1:end, :)=NaN;
handles.newtdats=newtdats;
set(handles.message_texts, 'String', 'all electrodes cleared')
update_ptcld_plot(hObject,handles)
% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% --- Executes on button press in fol_set_up.
function fol_set_up_Callback(hObject, eventdata, handles)
% hObject    handle to fol_set_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
seldir = uigetdir(cd, 'Select Directory with .ply files');
set(handles.message_texts,'String', 'Setting up Folders');
%get names of pointcloud
pcs=dir(seldir);
pcs={pcs(3:end).name};
subj_name=handles.subj_name;
mkdir([cd, '\dat\', subj_name, '\', subj_name, '_raw']);
dest=[cd, '\dat\', subj_name, '\', subj_name, '_raw'];
mkdir([cd, '\dat\', subj_name, '\', subj_name, '_saved'])
%move files
for j=1:length(pcs)
    source=[seldir, '\', pcs{j}];
    dest=[cd, '\dat\', subj_name, '\', subj_name, '_raw\', pcs{j}];
    copyfile(source, dest)
end
mkdir([cd, '\dat\', subj_name, '\', subj_name, '_elfound']);
mkdir([cd, '\dat\', subj_name, '\', subj_name, '_cropped']);
mkdir([cd, '\dat\', subj_name, '\', subj_name, '_pgdat']);
set(handles.message_texts,'String', 'Done Setting up Folders');
set(handles.fol_set_up, 'Visible', 'off')
guidata(hObject, handles);

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function fol_contrl_CreateFcn(hObject, ~, ~)
% hObject    handle to fol_contrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%--------------------------------------------------------------------------

function elthick_txt_Callback(hObject, eventdata, handles)
% hObject    handle to elthick_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elthick_txt as text
%        str2double(get(hObject,'String')) returns contents of elthick_txt as a double
subj_name=get(hObject,'String');
  uiresume;
%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function elthick_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elthick_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in folder_config.
function folder_config_Callback(hObject, eventdata, handles)
% hObject    handle to folder_config (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set background color
npx=200;
col=hsv2rgb([0, 0, .9])';
blank(:, :, 1)=ones(npx).*col(3);
blank(:, :, 2)=ones(npx).*col(3);
blank(:, :, 3)=ones(npx).*col(3);
imshow(blank)
handles.blank=blank;
guidata(hObject, handles)
set(handles.message_texts,'String', 'Type a subject name');
  set(handles.subj_name_btn, 'Visible', 'on')
uiwait
set(handles.message_texts,'String','Select folder with raw. ply files');
set(handles.subj_name_btn, 'Visible', 'off')
set(handles.fol_set_up, 'Visible', 'on')
uiwait;
guidata(hObject, handles)

%--------------------------------------------------------------------------
% --- Executes on button press in upload_ref_dat.
function upload_ref_dat_Callback(hObject, eventdata, handles)
% hObject    handle to upload_ref_dat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile;
%get subject name
if isfield(handles, 'subj_name')==1
subj_name=handles.subj_name;
else
    set(handles.message_texts, 'String', 'Type subject name')
    set(handles.subj_name_btn, 'Visible', 'on')
    uiwait;
    set(handles.subj_name_btn, 'Visible', 'off')
    subj_name=get(handles.subj_name_btn,'String');
end 
%set up pgdat folders if it doesnt already exist
folderName=[cd, '\dat\', subj_name, '\', subj_name, '_pgdat']
if isfolder(folderName)==0
%make directories 
mkdir([cd, '\dat\', subj_name, '\', subj_name, '_pgdat']);
end
%move reference file to the correct directory
base_pth=[cd, '\dat'];

    source=[path, file] ;
    dest=[base_pth,'\', subj_name, '\', subj_name '_pgdat'];
    %read and save as .mat file 
    file_ext=file(end-3:end);
switch file_ext
    case 'sfp'
pgdat = load_photogram_dat(path, file);
elpos=pgdat.elpos;
eval(['save ', [dest,'\coordinates.mat'], ' elpos']);
    case '..mat'
        copyfile(source, dest)
end
 set(handles.message_texts, 'String', 'Done setting up reference data')
guidata(hObject, handles);
%--------------------------------------------------------------------------
%BYE

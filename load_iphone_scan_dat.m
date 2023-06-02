%--------------------------------------------------------------------------
%
% Load all the iphone scanning data
%
%--------------------------------------------------------------------------
function [allscans,dat] = load_iphone_scan_dat(scnpth,scnfls,scnpth_raw,scnfls_raw)

%--------------------------------------------------------------------------
for n = 1:length(scnfls)
    %----------------------------------------------------------------------
    dat     = load([scnpth,'/',scnfls{n}]);
    if ~isempty(dat.colobj.Normal)
        disp('Loading Normals')
        %------------------------------------------------------------------
        elunlab = dat.newtdats;
        el = 90*pi/180;% str2double(get(handles.ini_rx,'String'))*pi/180;
        Rx          = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
        tdat        = (Rx*(dat.tdat'))';
        elunlab     = (Rx*(elunlab'))';        
        locspnvc    = (Rx*((dat.colobj.Location+dat.colobj.Normal)'))';    
        locs        = (Rx*(dat.colobj.Location'))';
        nvcs        = locspnvc - locs;                
        % De-mean
        cent        = mean(locs,1);
        tdat        = tdat    - repmat(cent,size(tdat,1),1);
        elunlab     = elunlab - repmat(cent,size(elunlab,1),1);
        locs        = locs    - repmat(cent,size(locs,1),1);
        %-----
        colobj      = pointCloud(locs, 'Color', dat.colobj.Color,'Normal',nvcs);
        
        %------------------------------------------------------------------
        allscans(n).colobj  = colobj;
        allscans(n).tdat    = tdat;
        allscans(n).elunlab = elunlab;        
        allscans(n).nvcs    = nvcs;
    else
        % disp('Loading Tris')
        %------------------------------------------------------------------
        % Load the triangles
        % 1628103088o778312oply_rawtris
        eval(['load ',scnpth_raw,'/',ifdec(scnfls{n}),'_rawtris', ...
            ' p t'])
        %------------------------------------------------------------------
        elunlab = dat.newtdats;
        el = 90*pi/180;% str2double(get(handles.ini_rx,'String'))*pi/180;
        Rx          = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
        tdat        = (Rx*(dat.tdat'))';
        elunlab     = (Rx*(elunlab'))';
        locs        = (Rx*(dat.colobj.Location'))';
        % De-mean
        cent        = mean(locs,1);
        tdat        = tdat    - repmat(cent,size(tdat,1),1);
        elunlab     = elunlab - repmat(cent,size(elunlab,1),1);
        locs        = locs    - repmat(cent,size(locs,1),1);
        %-----
        colobj      = pointCloud(locs, 'Color', dat.colobj.Color);
        
        %------------------------------------------------------------------
        allscans(n).colobj  = colobj;
        allscans(n).tdat    = tdat;
        allscans(n).elunlab = elunlab;
        allscans(n).tri     = t;
 
        figure
        trisurf(t,locs(:,1),locs(:,2),locs(:,3),'facecolor','cyan','linestyle','none')
        camlight left
    end
end
%-------------------------------------------------------------------------------
% 
% Construct a spherical mesh described by angles
% 
%-------------------------------------------------------------------------------
function angmsh = construct_sphang_msh(naz,nel,dbg_flg,topang)

if nargin < 4
    topang = 5;
elseif nargin < 3
    topang  = 5;
    dbg_flg = 0;
end
%-------------------------------------------------------------------------------
azs     = linspace(0,2*pi,naz);
els     = [-pi/2 linspace(-pi/2+topang*pi/180,pi/2-topang*pi/180,nel-2) pi/2];
angbds  = zeros((naz-1)*(nel-3)+2,4);
angcnts = zeros((naz-1)*(nel-3)+2,2);
k      = 1;
for n1 = 1:(naz-1)
    for n2 = 2:(nel-2)
        %-----------------------------------------------------------------------
        % Set the angle bounds
        angbds(k,:)  = [azs(n1) azs(n1+1)  els(n2) els(n2+1)];
        %-----------------------------------------------------------------------
        % Set the angle centers
        angcnts(k,:) = [mean(azs(n1:n1+1)) mean(els(n2:n2+1))];
        
        %-----------------------------------------------------------------------
        k = k+1;
    end
end
%-------------------------------------------------------------------------------
% Add the top and bottom caps
angbds(  k,:)  = [0 2*pi els(1)     els(2)];
angcnts(k,:)   = [0  els(1)];
angbds(k+1,:)  = [0 2*pi els(end-1) els(end)];
angcnts(k+1,:) = [0  els(end)];


%-------------------------------------------------------------------------------
% Construct a connection matrix 
L       = zeros((naz-1)*(nel-3)+2);
daz     = azs(2)-azs(1);
del     = els(2)-els(1);
for n = 1:(naz-1)*(nel-3)
    %---------------------------------------------------------------------------
    % the az+del point
    if angcnts(n,1) < max(angcnts(:,1))-1e-6     
        iazp = find( (abs( angcnts(:,1) - (angcnts(n,1)+daz))<1e-6) & (abs( angcnts(:,2) - angcnts(n,2))<1e-6) );
    else
        iazp = find( (abs( angcnts(:,1) - min(angcnts(:,1)))<1e-6) & (abs( angcnts(:,2) - angcnts(n,2))<1e-6) );
    end
    if ~isempty(iazp) 
        L(n,iazp) = -1;
        L(iazp,n) = -1;
        L(n,n)    = L(n,n)+1;
    end
    
    %---------------------------------------------------------------------------
    % the az-del point
    if angcnts(n,1) > min(angcnts(:,1))+1e-6     
        iazn = find( (abs( angcnts(:,1) - (angcnts(n,1)-daz))<1e-6) & (abs( angcnts(:,2) - angcnts(n,2))<1e-6) );
    else
        iazn = find( (abs( angcnts(:,1) - max(angcnts(:,1)))<1e-6) & (abs( angcnts(:,2) - angcnts(n,2))<1e-6) );
    end    
    if ~isempty(iazp) 
        L(n,iazn) = -1;
        L(iazn,n) = -1;
        L(n,n)    = L(n,n)+1;
    end
    %---------------------------------------------------------------------------
    % the el+del point
    if angcnts(n,2) < max(angcnts(1:end-2,2))-1e-6     
        ielp = find( (abs( angcnts(:,2) - (angcnts(n,2)+del))<1e-6) & (abs( angcnts(:,2) - angcnts(n,2))<1e-6) );
    else
        ielp = find( (abs( angcnts(:,2) - max(angcnts(:,2)))<1e-6) );
    end
    if ~isempty(ielp) 
        L(n,ielp) = -1;
        L(ielp,n) = -1;
        L(n,n)    = L(n,n)+1;
    end
    
    %---------------------------------------------------------------------------
    % the el-del point
    if angcnts(n,2) > min(angcnts(1:end-2,2))+1e-6     
        ieln = find( (abs( angcnts(:,2) - (angcnts(n,2)-del))<1e-6) & (abs( angcnts(:,1) - angcnts(n,1))<1e-6) );
    else
        ieln = find( (abs( angcnts(:,2) - min(angcnts(:,2)))<1e-6) );
    end
    if ~isempty(iazp) 
        L(n,ieln) = -1;
        L(ieln,n) = -1;
        L(n,n)    = L(n,n)+1;
    end
        
end



%-------------------------------------------------------------------------------
angmsh.angbds  = angbds;
angmsh.angcnts = angcnts;
angmsh.L       = L;


%-------------------------------------------------------------------------------
% Plot the angular mesh
if dbg_flg == 1
    plot_angmesh(angmsh)
end
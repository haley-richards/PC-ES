%--------------------------------------------------------------------------
%
% Compute the normal of each electrode
%
%--------------------------------------------------------------------------
function allscans = comp_elec_nvecs(allscans,elthick,dbg_flg)

%--------------------------------------------------------------------------
% Set parameters
elrad = 3/1000; % radius

%--------------------------------------------------------------------------
% Loop through each scan, make the normal vector, align them, and
av_num_tris = zeros(length(allscans),1);
for n = 1:length(allscans)
    %----------------------------------------------------------------------
    tdat    = allscans(n).tdat;
    nv_tdat = NaN*ones(size(tdat));
    if isfield(allscans(n),'tri')
        p       =double(allscans(n).colobj.Location);
        t       = allscans(n).tri;
        %p=p*1000;
        tcs     = double(get_tcs(p,t));
        %tcs=tcs*1000;
        %----------------------------------------------------------------------
        p=double(p);
        a = p(t(:,2),:) - p(t(:,1),:);
        b = p(t(:,3),:) - p(t(:,1),:);
        % the normal is the cross product of two vectors
        % a x b = (a2*b3-a3*b2)i + (a3*b1-a1*b3)j + (a1*b2-a2*b1)k
        ns = [ ...
            (a(:,2).*b(:,3) - a(:,3).*b(:,2)) ...
            (a(:,3).*b(:,1) - a(:,1).*b(:,3)) ...
            (a(:,1).*b(:,2) - a(:,2).*b(:,1))];
 
        % Normalize
        nlens = sqrt(sum(ns.^2,2));
        nlens(nlens==0)=.00000001;
        ns(ns==0)=.000000001;
        ns = ns ./ repmat(nlens,1,3);
        
        %------------------------------------------------------------------
        % Make sure the normal vector point away from the normal vector
        % dot the triangle center with the normal vector, if they are negative,
        % then they are not aligned correctly
        dots   = sum(tcs.*ns,2);
        iflips = find( dots < 0);
        ns(iflips,:) = -ns(iflips,:);
        
    elseif isfield(allscans(n),'nvcs')
        %------------------------------------------------------------------
        % If there already are normal vectors, then just use them.
        %------------------------------------------------------------------
        tcs = allscans(n).colobj.Location;
        ns  = allscans(n).nvcs;
        nlens = sqrt(sum(ns.^2,2));
        ns = ns ./ repmat(nlens,1,3);
        
    end
    %------------------------------------------------------------------
    % Loop through the tdat and construct the normal vectors
    numtris = [];
    for k = 1:size(tdat,1)
        if isnan(tdat(k,1)) == 0
            is           = find( sqrt(sum( (tcs - tdat(k,:)).^2,2)) < elrad);
            nv_tdat(k,:) = mean( ns(is,:),1);
            nv_tdat(k,:) = nv_tdat(k,:) / norm(nv_tdat(k,:));
            numtris      = [numtris; length(is)];
        end
    end
    av_num_tris(n) = mean(numtris);
    
    
    if dbg_flg == 1
        disp(['Ave tri to make normal: ',num2str(av_num_tris(n))])
        % Check the normal vectors
        MS = 200; % Marksize of point cloud
        figure
        pcshow(allscans(n).colobj,'Markersize',MS);
        hold on
        %------------------------------------------------------------------
        % Labeled electrodes
        plot3(tdat(:,1),tdat(:,2),tdat(:,3),'.y','markersize',36)
        fac = 10/1000;
        for k = 1:size(tdat,1)
            if isnan(tdat(k,1)) == 0
                plot3(tdat(k,1) + fac*[0 nv_tdat(k,1)], ...
                    tdat(k,2) + fac*[0 nv_tdat(k,2)], ...
                    tdat(k,3) + fac*[0 nv_tdat(k,3)],'-r')
            end
        end
        %--------------------------------------------------------------------------
        % Formatting
        axis equal
        grid on
        box on
        FS = 12;
        xlabel('X','fontsize',FS,'FontName','times')
        ylabel('Y','fontsize',FS,'FontName','times')
        zlabel('Z','fontsize',FS,'FontName','times')
        set(gca,'FontSize',FS,'FontName','times')
        
    end
    
    %----------------------------------------------------------------------
    % Move the electrode back by its thickness
    for k = 1:size(tdat,1)
        if isnan(tdat(k,1)) == 0
            tdat(k,:) = tdat(k,:) - elthick*nv_tdat(k,:);
            
        end
    end
    
    %----------------------------------------------------------------------
    % Record the normal vectors and the updated tdat points
    allscans(n).nv_tdat = nv_tdat;
    allscans(n).tdat    = tdat;
    
end


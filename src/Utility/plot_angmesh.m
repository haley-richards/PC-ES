%-------------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------------
function plot_angmesh(angmsh,is_reg,plt_typ,rad)

if nargin == 1
    is_reg  = 1:size(angmsh.angbds,1);
    plt_typ = [1 1];
    rad     = 1;
end

%--------------------------------------------------------------------------
angbds = angmsh.angbds;
Ns     = 30;

if plt_typ(1) == 1
    % figure
    hold on
    pcols = plot_cols();
    k     = 1;
    for n = 1:size(angbds,1)
        if ~isempty(intersect(is_reg,n))
            %--------------------------------------------------------------
            % Sample the boundary
            azs = linspace(angbds(n,1),angbds(n,2),Ns);
            els = linspace(angbds(n,3),angbds(n,4),Ns);
            %--------------------------------------------------------------
            fill3( ...
                rad*[cos(azs)*cos(els(1)) cos(azs(end))*cos(els) cos(azs(end:-1:1))*cos(els(end)) cos(azs(1))*cos(els(end:-1:1))], ...
                rad*[sin(azs)*cos(els(1)) sin(azs(end))*cos(els) sin(azs(end:-1:1))*cos(els(end)) sin(azs(1))*cos(els(end:-1:1))], ...
                rad*[sin(els(1))+0*azs     sin(els)               sin(els(end))+0*azs              sin(els(end:-1:1))], ...
                pcols{k},'facealpha',0.5)
            k = k+1;
            if k > length(pcols)
                k = 1;
            end
        end
    end
    
    %----------------------------------------------------------------------
    % Plot the element centers
    angcnts = angmsh.angcnts;
    cnts    = rad*[cos(angcnts(:,1)).*cos(angcnts(:,2)) sin(angcnts(:,1)).*cos(angcnts(:,2)) sin(angcnts(:,2))];
    plot3(cnts(is_reg,1),cnts(is_reg,2),cnts(is_reg,3),'.r','markersize',22)
    view(3)
    title(['Num. elements: ',num2str(size(angcnts,1))])
end

%--------------------------------------------------------------------------
% Plot the connections
if plt_typ(2) == 1
    
    L = angmsh.L;
    for n = 1:size(L,1)
        is = find( L(n,:) < 0);
        for k = 1:length(is)
            s = 0.3;
            plot3([cnts(n,1) (1-s)*cnts(n,1) + s*cnts(is(k),1)], ...
                [cnts(n,2) (1-s)*cnts(n,2) + s*cnts(is(k),2)], ...
                [cnts(n,3) (1-s)*cnts(n,3) + s*cnts(is(k),3)],'-b')
        end
    end
    figure
    imagesc(L)
    
end


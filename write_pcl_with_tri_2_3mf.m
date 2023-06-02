%-------------------------------------------------------------------------------
%
% Write a surface triangulation and a vector to be displayed as a bunch of
% rectangle patch objects to an STL file
%
%-------------------------------------------------------------------------------
function write_pcl_with_tri_2_3mf(objsavnme,plt_ts,plt_ps,plt_cols, ...
    rectsxs,rectsys,rectszs,svecs, ...
    clim,cmap,dots,dot_cs,dot_diam)

%-------------------------------------------------------------------------------
% Convert the patch object into a surface triangulation
c_ptch_tot = [];
t_ptch_tot = [];
p_ptch_tot = [];
for k = 1:length(rectsxs)
    rectsx = rectsxs{k};
    rectsy = rectsys{k};
    rectsz = rectszs{k};
    svec   = svecs{k};
    t_ptch = zeros(2*size(rectsx,2),3);
    p_ptch = zeros(4*size(rectsx,2),3);
    c_ptch = zeros(2*size(rectsx,2),3);
    for n = 1:size(rectsx,2)
        p_ptch(4*(n-1)+(1:4),1) = rectsx(:,n);
        p_ptch(4*(n-1)+(1:4),2) = rectsy(:,n);
        p_ptch(4*(n-1)+(1:4),3) = rectsz(:,n);
        t_ptch(2*(n-1)+1,:) = 4*(n-1) + [1 2 3];
        t_ptch(2*(n-1)+2,:) = 4*(n-1) + [1 3 4];
        ind  = round(size(cmap,1)/(clim(2)-clim(1))*(svec(n)-clim(1)));
        if ind > size(cmap,1)
            ind = size(cmap,1);
        elseif ind < 1
            ind = 1;
        end
        c_ptch(2*(n-1)+1,:) = cmap(ind,:);
        c_ptch(2*(n-1)+2,:) = cmap(ind,:);
    end
    [t_ptch,p_ptch] = remove_repeated_nodes(t_ptch,p_ptch);
    
    % Combine the triangulations
    c_ptch_tot = [c_ptch_tot; c_ptch];
    t_ptch_tot = [t_ptch_tot; t_ptch+size(p_ptch_tot,1)];
    p_ptch_tot = [p_ptch_tot; p_ptch];
end
p_ptch = p_ptch_tot;
t_ptch = t_ptch_tot;
c_ptch = c_ptch_tot;

%-------------------------------------------------------------------------------
% Construct a a surface triangulation for all the dots
[X,Y,Z] = sphere(4);
sph_p   = dot_diam/2*[X(:) Y(:) Z(:)];
sph_t   = convhulln(sph_p);
dot_p = [];
dot_t = [];
dot_c = [];
for n = 1:size(dots,1)
    dot_t = [dot_t; sph_t+size(dot_p,1)];
    dot_p = [dot_p; sph_p + repmat(dots(n,:),size(sph_p,1),1)];
    dot_c = [dot_c; repmat(dot_cs(n,:),size(sph_t,1),1)];
end
% whos
% size(dots,1)
% return
% clear

%-------------------------------------------------------------------------------
% Combine the plt_ts, plt_ps triangulations
p = [];
t = [];
cmat = [];
for n = 1:length(plt_ts)
    t = [t; plt_ts{n}+size(p,1)];
    p = [p; plt_ps{n}];
    cmat = [cmat;  repmat(plt_cols{n},size(plt_ts{n},1),1)];
end

%-------------------------------------------------------------------------------
% Combine the triangulations
% cmat  = [ repmat([0 1 1],size(t,1),1); c_ptch; dot_c];
% cmat  = [ repmat(0.5*[1 1 1],size(t,1),1); c_ptch; dot_c];
cmat  = [ cmat; c_ptch; dot_c];
t     = [t; t_ptch+size(p,1); dot_t+size(p,1)+size(p_ptch,1)];
p     = [p; p_ptch; dot_p];
[t,p] = remove_unused_nodes(t,p);
% Convert the color to each vertex
DT     = triangulation(t,p);
iss    = vertexAttachments(DT);
cmatv  = zeros(size(p,1),3);
for k = 1:size(p,1)
    cmatv(k,:) = mean(cmat(iss{k},:),1);
end



%-------------------------------------------------------------------------------
% figure
% hold on
% trisurf(t,p(:,1),p(:,2),p(:,3),sqrt(sum(cmat.^2,2)),'facealpha',0.5)
%-------------------------------------------------------------------------------
% WRite the stl file
% fv.vertices = p;
% fv.faces    = t;
% fv
% stlwrite('testfile.stl',fv,'MODE','binary', ...
%     'FACECOLOR',round(255*cmat+1));%     - Single colour (1-by-3) or one-colour-per-face (N-by-3) 
addpath S:\digihisto\Ethan\software\cvergari-write3mf-1136acb
 write3mf([objsavnme,'.3mf'], p , t, uint8(round(255*cmatv)))
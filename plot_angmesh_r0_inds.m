%-------------------------------------------------------------------------------
% 
% 
% 
%-------------------------------------------------------------------------------
function plot_angmesh_r0_inds(angmsh,r0,inds)

%-------------------------------------------------------------------------------
angbds = angmsh.angbds;
Ns     = 30;

for i = 1:length(inds)
    n  = inds(i);
    %---------------------------------------------------------------------------
    % Sample the boundary     
    azs = linspace(angbds(n,1),angbds(n,2),Ns);
    els = linspace(angbds(n,3),angbds(n,4),Ns);
    %---------------------------------------------------------------------------
    fill3( ...
        r0*[cos(azs)*cos(els(1)) cos(azs(end))*cos(els) cos(azs(end:-1:1))*cos(els(end)) cos(azs(1))*cos(els(end:-1:1))], ...
        r0*[sin(azs)*cos(els(1)) sin(azs(end))*cos(els) sin(azs(end:-1:1))*cos(els(end)) sin(azs(1))*cos(els(end:-1:1))], ...
        r0*[sin(els(1))+0*azs     sin(els)               sin(els(end))+0*azs              sin(els(end:-1:1))], ...
        'r','facealpha',0.25)     
end

%-------------------------------------------------------------------------------
% 
% Calculate the area of the sectors (not including a radius factor (r^2).
% So its either assumed to be constant and will be incorporated later, or
% this is if r = 1)
% 
%-------------------------------------------------------------------------------
function secas = calc_angmsh_areas(angmsh)

angbds = angmsh.angbds;
nsec   = size(angbds,1);
secas  = zeros(nsec,1);
for n = 1:nsec
    azb = sort(angbds(n,1:2));
    phb = sort(pi/2-angbds(n,3:4)); % defined 0 at +z
    secas(n) = (azb(2)-azb(1))*-(cos(phb(2))-cos(phb(1)));
end

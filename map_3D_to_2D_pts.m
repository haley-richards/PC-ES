%--------------------------------------------------------------------------
%
% Convert from 3D points to a 2D plane. The transformation relies on a
% center (Porg) and a horizontal plane at a fixed z-value. Each point is
% connected to the plane via a line: 
%   r(t) = [(x-Porg(1))*t+Porg(1),(y-Porg(2))*t+Porg(2),(z-Porg(3))*t+Porg(3)]
% 
% where t is such that r(t)_z = the fixed z-value. So, t is given by
%   t = (zpln-Porg(3))/(z-Porg(3));
%  
% How do we go back?
%   Given the projection points x_p,y_p we want to know x,y,z
%   x_p = (x-Porg(1))*t + Porg(1)
%   y_p = (y-Porg(2))*t + Porg(2)
%   or
%   x = (x_p - Porg(1))/t + Porg(1) 
%   y = (y_p - Porg(2))/t + Porg(2)
%   z = (zpln- Porg(3))/t + Porg(3) 
% 
%--------------------------------------------------------------------------
function elps_2d = map_3D_to_2D_pts(elps,Porg)


zpln     = 1;
rays     = elps - repmat(Porg,size(elps,1),1);
t        = (zpln - Porg(3))./rays(:,3); 
elps_2d  = [rays(:,1).*t + Porg(1) rays(:,2).*t + Porg(2)];

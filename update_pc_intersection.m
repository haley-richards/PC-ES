function [pcB_out]=update_pc_intersection(hobject, handles, n_1, n_2) 
% get sub indices of first point cloud (that has been updated)
switch n_1
    case 1
        A=handles.pc_sub1;
        pcA=handles.pc1;
    case 2
        A=handles.pc_sub2;
         pcA=handles.pc2;
    case 3
        A=handles.pc_sub3;
         pcA=handles.pc3;
    case 4
        A=handles.pc_sub4;
         pcA=handles.pc4;
end

% get sub indices of second point cloud (that needs to be updates)
switch n_2
    case 1
        B=handles.pc_sub1;
         pcB=handles.pc1;
    case 2
        B=handles.pc_sub2;
         pcB=handles.pc2;
    case 3
        B=handles.pc_sub3;
         pcB=handles.pc3;
    case 4
        B=handles.pc_sub4;
         pcB=handles.pc4;
end
%locsA=pcA.Location;
colsA=pcA.Color;

locsB=pcB.Location;
colsB=pcB.Color;
%update the intersersecting regions of subscans
[int, ina, inb]=intersect(A, B);
if isempty(int)==0
int_cols_a=colsA(ina, :);
%int_cols_b=colsB(inb, :)

%update colors of intersecting pointcloud
colsB(inb, :)=int_cols_a;
end
%create new point cloud
pcB_out=pointCloud(locsB, 'Color', colsB);

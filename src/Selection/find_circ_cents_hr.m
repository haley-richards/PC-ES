%-------------------------------------------------------------------------------
%
% Run circle center finding algorithm
%
%-------------------------------------------------------------------------------
function [centers] = find_circ_cents_hr(img_m,pts_m,dbg_flg)
% save testdata
% clear

% clear
% clc
% close all
% load testdata

%---------------------------------------------------------------------------
% Plot the image and pick an electrode as the template
figure
set(gcf,'position',[257         124        1198         825])
imagesc(img_m)
hold on
if size(pts_m,1) > 0
    plot(pts_m(:,1),pts_m(:,2),'.r','markersize',12)
end
colormap gray
hold on
axis equal
set(gca, 'YDir','normal')
title('Pick unmarked electrodes, when done hit the ESC (escape) button')
% butt    = 1;
% centers = [];
% while (butt == 1) || (butt == 2)
%     [x1,y1,butt] = ginputWhite(1);
%     if butt == 1        
%         plot(x1,y1,'.g','markersize',18)
%         centers = [centers; x1,y1];
%     elseif butt == 2
%         [tmp,i] = min( sum( (centers - repmat([x1 y1],size(centers,1),1)).^2,2));
%         centers(i,:) = [];
%         clf
%         imagesc(img_m)
%         hold on
%         if size(pts_m,1) > 0
%             plot(pts_m(:,1),pts_m(:,2),'.r','markersize',18)
%         end
%         if size(centers,1) > 0 
%             plot(centers(:,1),centers(:,2),'.g','markersize',18)
%         end        
%         % colormap gray
%         hold on
%         axis equal
%         set(gca, 'YDir','normal')        
%     end
% end
select_flag=1;
pt_count=1;
while select_flag==1
select=drawpoint;
if isempty(select.Position)==0
centers(pt_count, :)=select.Position;
%disp([num2str(pt_count),' point(s) recorded in this sub scan!'])
pt_count=pt_count+1;
plot(select.Position(1), select.Position(2), 'r.', 'MarkerSize', 30) 
else
    select_flag=0;
end
end



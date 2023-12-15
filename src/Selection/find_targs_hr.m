%-------------------------------------------------------------------------------
%
% Run circle center finding algorithm
%
%-------------------------------------------------------------------------------
function [centers, targ_label] = find_targs_hr(img_m,dbg_flg)
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
% hold on
% if size(pts_m,1) > 0
%     plot(pts_m(:,1),pts_m(:,2),'.r','markersize',12)
% end
colormap gray
hold on
axis equal
set(gca, 'YDir','normal')
title('Pick unmarked electrodes, Hit ESC (escape) button when done.')
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
label_cnt=1;
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
%GET TARGET POINT COLORS
%create matlab listbox
if pt_count>1
    fig=uifigure;
    label1 = uilabel(fig,...
        'Position',[50 150 500 200],...
        'Text','type target number then click outside of box (in order clicked).');
    label2 = uilabel(fig,...
        'Position',[50 100 300 150],...
        'Text', 'Use nominal map for reference (electrodes in red).');
    textarea = uitextarea(fig,...
        'Position',[100 100 150 60],...
        'ValueChangedFcn',@(textarea,event) textEntered(textarea, label2));
    %-----------------------------------------------------------------------------
    %get user data from listbox
    for i=1:size(centers, 1);
        uiwait; %wait until user has entered text in the window
        list_label{i}=textarea.Value;
        targ_label(label_cnt)=[str2num(cell2mat(list_label{i}))];
        label_cnt=label_cnt+1;
        %------------------------------------------------------------------------------------
    end
else
    centers=[];
    targ_label=[];
end
hold off
fig.Visible='off'
end


%------------------------------------------------------------------------
% listbox entry function
function textEntered(textarea,label2, label1)
val = textarea.Value;
% Check each element of text area cell array for text
for k = 1:length(val)
    if(~isempty(val{k}))
        label1.Text = 'labeled! Go to the next scan';
        label2.Text=['You entered ', cell2mat(val)];
        uiresume
        break;
    end
end
end



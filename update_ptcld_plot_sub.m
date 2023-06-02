function  update_ptcld_plot_sub(hObject,handles)
MS     = str2double(get(handles.MStag,'String'));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Update the figure
colobj=handles.subscn;
axes(handles.mainax);
azel = get(handles.mainax,'view');
% azel=handles.proj_data;
% azel=double(azel(handles.scn_num, 1:2));
% tic
cla
% disp(['PC update toc 1:',num2str(toc)])
% tic
pcshow(colobj,'Markersize',MS);
% disp(['PC update toc 2:',num2str(toc)])
hold on
if isfield(handles,'tdat') == 1
%    % plot3(tdat(:,1),tdat(:,2),tdat(:,3),'.y','markersize',36)
%     popstr = get(handles.popupmenu1,'String');
%     
%     for n = 1:length(popstr)
%         if (length(popstr{n}) < 4) && (isnan(tdat(n,1)) == 0)
%             tmp = double(tdat(n,:));
%             text(1.02*tmp(1),1.02*tmp(2),1.02*tmp(3),popstr{n},'fontsize',16,'FontName','times','color','red')
%         end
%     end
end
% if isfield(handles,'eclks') == 1
%     plot3(eclks(:,1),eclks(:,2),eclks(:,3),'.m','markersize',32)
% end

% % if handles.plotplane == 1
% %     xbs = [str2double(get(handles.zpl_x1,'String')) ...
% %         str2double(get(handles.zpl_x2,'String'))];
% %     ybs = [str2double(get(handles.zpl_y1,'String')) ...
% %         str2double(get(handles.zpl_y2,'String'))];
% %     z0  = str2double(get(handles.zpl_z,'String'));
% %     fill3([xbs(1) xbs(2) xbs(2) xbs(1)],[ybs(1) ybs(1) ybs(2) ybs(2)],z0*[1 1 1 1],'red')
% % end
% % camlight headlight
%--------------------------------------------------------------------------
% tic
% if isfield(handles,'add_bestfit_els') == 1
%     if handles.add_bestfit_els == 1
%         set(handles.mainax,'view',azel);
%         [Bout,elns,bestsclf] = try_bestfitting_of_nominal_elecs(handles);
%         handles.Bout = Bout;
%         handles.elns = elns;        
%         handles.bestsclf = bestsclf;
%     end
% end
% disp(['PC update toc 3:',num2str(toc)])
%--------------------------------------------------------------------------
% tic
axis equal
grid on
box on
FS = 12;
xlabel('X','fontsize',FS,'FontName','times')
ylabel('Y','fontsize',FS,'FontName','times')
zlabel('Z','fontsize',FS,'FontName','times')
set(gca,'FontSize',FS,'FontName','times')
set(handles.mainax,'view',azel);
% disp(['PC update toc 4:',num2str(toc)])
%--------------------------------------------------------------------------
% Update handles structure
% tic
guidata(hObject, handles);
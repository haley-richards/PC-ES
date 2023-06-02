function  update_small_axis(hObject,handles)
MS     = str2double(get(handles.MStag,'String'));
%--------------------------------------------------------------------------
% Perform translations of the points
% if handles.cropped == 1
    colobj=handles.colobj_orig;
% else
%     colobj = handles.colobj;
%     colobj = pcdownsample(colobj,'random',0.01);
%     set(handles.message_texts,'String','Please Crop before rotating/translating');
% end
%--------------------------------------------------------------------------
% Update the figure
axes(handles.axes1);
if isfield(handles, 'is_sub')==1
if handles.is_sub ==1
       scn_num=handles.scn_num;
            switch scn_num
                case 1
                    col_inds=handles.pc_sub1;
              case 2
                     col_inds=handles.pc_sub2;
              case 3
                     col_inds=handles.pc_sub3;
              case 4
                     col_inds=handles.pc_sub4;
            end
   else
    pcshow(handles.colobj);         
end

end
cols=colobj.Color;
locs=colobj.Location;
cols(col_inds, :)=repmat([0, 250, 0], length(col_inds), 1);
colobj_view=pointCloud(locs,'Color', cols);
pcshow(colobj_view)
% azel = get(handles.mainax,'view');
% % tic
% cla
% disp(['PC update toc 1:',num2str(toc)])
% tic
% im=imread('EGI_black_electrode_map.jpg');
% imshow(im)
% title('reference map (black electrodes in red)')
%pcshow(colobj,'Markersize',MS);
% disp(['PC update toc 2:',num2str(toc)])
%hold on
% if isfield(handles,'tdat') == 1
%     plot3(tdat(:,1),tdat(:,2),tdat(:,3),'.y','markersize',36)
% end
%     popstr = get(handles.popupmenu1,'String');
%     
%     for n = 1:length(popstr)
%         if (length(popstr{n}) < 4) && (isnan(tdat(n,1)) == 0)
%             tmp = double(tdat(n,:));
%             text(1.02*tmp(1),1.02*tmp(2),1.02*tmp(3),popstr{n},'fontsize',16,'FontName','times','color','red')
%         end
%     end
% % end
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
%set(handles.mainax,'view',azel);
% disp(['PC update toc 4:',num2str(toc)])
%--------------------------------------------------------------------------
% Update handles structure
% tic
axes(handles.mainax);
guidata(hObject, handles);
% disp(['PC update toc 5:',num2str(toc)])
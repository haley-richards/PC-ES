%-------------------------------------------------------------------------------
%
% Calculate the overall error
%
%-------------------------------------------------------------------------------
function valtrans = calc_valerr(pgdat,tdat_mrg,colobj_mrg,dbg_flg)


is_fnd = find(isnan(tdat_mrg(:,1))==0);

disp([num2str(length(is_fnd)),' Electrodes found'])

%-------------------------------------------------------------------------------
% Collect the target and corresponding nominal electrode points
A  = [tdat_mrg(is_fnd,:) ones(length(is_fnd),1)]';
B  = [pgdat.elpos(is_fnd,:)/100 ones(length(is_fnd),1)]';

figure
hold on
plot3(tdat_mrg(is_fnd,1),tdat_mrg(is_fnd,2),tdat_mrg(is_fnd,3),'.k','markersize',16)
plot3(pgdat.elpos(is_fnd,1)/100,pgdat.elpos(is_fnd,2)/100,pgdat.elpos(is_fnd,3)/100,'.r','markersize',16)

%-------------------------------------------------------------------------------
% Get registered points from each frame: T is such that A = T * B
[T, rmserr]   = transform_loc(A, B, 1);
disp(['RMS error: ',num2str(rmserr*1000),' mm'])
T
valtrans = (T*([pgdat.elpos/100 ones(size(pgdat.elpos,1),1)]'))';

figure
hold on

%-------------------------------------------------------------------------------
if dbg_flg == 1
        
    %---------------------------------------------------------------------------
    % Plot the nominal eeg cap and the 
    figure
    set(gcf,'position',[680         394        1109         584])
    subplot(1,2,1)
    hold on
    plot3(tdat_mrg(:,1),tdat_mrg(:,2),tdat_mrg(:,3),'.k','markersize',16)
    plot3(valtrans(:,1),valtrans(:,2),valtrans(:,3),'.r','markersize',16)   
    for n = 1:size(tdat_mrg,1)    
        plot3([tdat_mrg(n,1) valtrans(n,1)],[tdat_mrg(n,2) valtrans(n,2)],[tdat_mrg(n,3) valtrans(n,3)],'-g')
    end
    axis equal
    grid on
    box on
    FS = 12;
    xlabel('X','fontsize',FS,'FontName','times')
    ylabel('Y','fontsize',FS,'FontName','times')
    zlabel('Z','fontsize',FS,'FontName','times')
    set(gca,'FontSize',FS,'FontName','times')
    set(gcf,'color','w');
    legend('iPhone','Photogrammetry')
    axis off
    
    subplot(1,2,2) % figure
    hold on
    MS = 100;        
    pcshow(colobj_mrg,'Markersize',MS);
    axis equal
    grid on
    box on
    FS = 12;
    xlabel('X','fontsize',FS,'FontName','times')
    ylabel('Y','fontsize',FS,'FontName','times')
    zlabel('Z','fontsize',FS,'FontName','times')
    set(gca,'FontSize',FS,'FontName','times')
    set(gcf,'color','w');
    axis off
    
    %     azs = [0:45:(360-45)]
    %     els = [15 30];
    %     for n2 = 1:length(els)
    %         for n1 = 1:length(azs)
    %             set(gca,'view',[azs(n1) els(n2)]);
    %
    %         end
    %     end
    
    
    writerObj = VideoWriter('figs/mergefinal_pointcloud_mov.mp4','MPEG-4');
    writerObj.FrameRate = 2;
    open(writerObj);
    
    azs   = (0:5:360)';
    azels = [ ...
        azs 22.5+0*(0:5:360)'; ...
        azs 45+0*(0:5:360)'; ...
        ];
    for k = 1:size(azels,1)
        for j = 1:2
            subplot(1,2,j)
            view(azels(k,:))
        end
        %-------------------------------------------------------------------
        F = getframe(gcf);
        writeVideo(writerObj,F);
        pause(0.01)
        
    end
    close(writerObj);
    
    
end









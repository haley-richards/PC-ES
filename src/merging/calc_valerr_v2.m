%--------------------------------------------------------------------------
%
% Calculate the overall error between the photogrammetry data and the
% nominal cap data. 
%
%--------------------------------------------------------------------------
function nomps0 = calc_valerr_v2(nomdat,tdat_mrg,colobj_mrg,dbg_flg,savnme)


is_fnd = find(isnan(tdat_mrg(:,1))==0);

disp([num2str(length(is_fnd)),' Electrodes found'])

%--------------------------------------------------------------------------
% 1. Perform a best-fit between the nominal and registration data.
nomps0 = bestfit_nomel_noplot_fullset_loc(tdat_mrg,nomdat);

%--------------------------------------------------------------------------
% Get registered points from each frame: T is such that A = T * B
N = min([size(tdat_mrg,1) size(nomps0,1)]);
ds_nom2tdat = sqrt(sum( (tdat_mrg(1:N,:) - nomps0(1:N,:)).^2,2));
rmserr      = sqrt(1/length(is_fnd)*sum(ds_nom2tdat(is_fnd).^2 ));

disp(['RMS error: ',num2str(rmserr*1000),' mm'])
disp(['   Max diff: ',num2str(round(max(ds_nom2tdat,[],'omitnan')*1000,1)),' mm'])
disp(['  Mean diff: ',num2str(round(mean(ds_nom2tdat,'omitnan')*1000,1)),' mm'])
disp(['Median diff: ',num2str(round(median(ds_nom2tdat,'omitnan')*1000,1)),' mm'])

figure;set(gcf,'position',[297         348        1243         483])
subplot(1,2,1); hold on
plot(ds_nom2tdat*1000,'.k','markersize',16)
lbl_fmt_fig('Electrode','Distance (mm)','Distance to iPhone Scan','','',12)
legend('Nominal','Photogrammetry')
subplot(1,2,2);
histogram(ds_nom2tdat(is_fnd)*1000)
lbl_fmt_fig('Distance (mm)','Frequency','Nominal vs iPhone Scan','','',12)

%--------------------------------------------------------------------------
if dbg_flg == 1
        
    %---------------------------------------------------------------------------
    % Plot the nominal eeg cap and the 
    figure
    set(gcf,'position',[680         394        1109         584])
    subplot(1,2,1)
    hold on
    plot3(tdat_mrg(:,1),tdat_mrg(:,2),tdat_mrg(:,3),'.k','markersize',16)
    plot3(nomps0(:,1),nomps0(:,2),nomps0(:,3),'.r','markersize',16)   
    for n = 1:size(tdat_mrg,1)    
        plot3([tdat_mrg(n,1) nomps0(n,1)],[tdat_mrg(n,2) nomps0(n,2)],[tdat_mrg(n,3) nomps0(n,3)],'-g')
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
    legend('iPhone','Nominal')
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
    
    
    writerObj = VideoWriter([savnme,'_mergefinal_pointcloud_mov.mp4'],'MPEG-4');
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









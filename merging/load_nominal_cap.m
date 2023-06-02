%-------------------------------------------------------------------------------
%
% Load the nominal cap
%
%-------------------------------------------------------------------------------
function nomdat = load_nominal_cap(dbg_flg)

load dat/UCL_head_mesh_256eeg_elecs msh el256 eeg_mps
nomdat.msh     = msh;
nomdat.el256   = el256;
nomdat.eeg_mps = eeg_mps/100; % convert from cm to mm
fnm     = 'netModel_ideal.xml';
mdl     = read_docbrown_file(['dat\',fnm]);

%---------------------------------------------------------------------------
% Collect the target and corresponding nominal electrode points
A  = [nomdat.eeg_mps ones(size(nomdat.eeg_mps,1),1)]';
B  = [mdl.sens_ps/100 ones(258,1)]';
[T, rmserr]   = transform(A, B(:,1:256), 1);
eeg_mps2 = (T*B)';

nomdat.eeg_mps  = [nomdat.eeg_mps; ...
    mean(nomdat.eeg_mps([9 45 81 132 186]',:),1); ...
    mean(nomdat.eeg_mps([89 90 130 129 101 100]',:),1)];
nomdat.eeg_mps2 = eeg_mps2(:,1:3);
nomdat.cnxs     = mdl.cnxs;

if dbg_flg == 1
    figure
    hold on
    plot3(nomdat.eeg_mps(:,1),nomdat.eeg_mps(:,2),nomdat.eeg_mps(:,3),'.k','markersize',12)
    for n = 1:size(nomdat.eeg_mps,1)
        text(nomdat.eeg_mps(n,1),nomdat.eeg_mps(n,2),nomdat.eeg_mps(n,3),num2str(n),'fontsize',16)
    end
    plot3(nomdat.eeg_mps(256,1),nomdat.eeg_mps(256,2),nomdat.eeg_mps(256,3),'.b','markersize',22)
    plot3(nomdat.eeg_mps(257,1),nomdat.eeg_mps(257,2),nomdat.eeg_mps(257,3),'.g','markersize',22)
    
    figure
    hold on
    plot3(nomdat.eeg_mps(:,1),nomdat.eeg_mps(:,2),nomdat.eeg_mps(:,3),'.k','markersize',12)
    plot3(eeg_mps2(:,1),eeg_mps2(:,2),eeg_mps2(:,3),'.r','markersize',12)
    legend('Brainstorm','photogrammetry system')
    
    figure
    hold on
    mdl
    mdl.sens_ps = nomdat.eeg_mps*100;
%     mdl.sens_ps = nomdat.eeg_mps2*100;
    plot_docbrown(mdl,'b')
end

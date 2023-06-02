%-------------------------------------------------------------------------------
% 
% merge scans and autolabel
% 
%-------------------------------------------------------------------------------
function [rms_err, rms_fullset]= comp_pg_tdat_vf(tdat_mrg, colobj_mrg, valpth, subj_name, file_ext, savepth)

% Load the nominal cap
nomdat = load_nominal_cap(0);
%-------------------------------------------------------------------------------
%load photogrammetry data
% Load the photogrammetry data
if file_ext=='.sfp'
pgdat = load_photogram_dat(valpth,'coordinates.sfp');
else
    pgdat = load_matlab_pgdat(valpth);
end
%-------------------------------------------------------------------------------
%compare merged coordinates with photogrammetry 
shrt_name=subj_name;
[valtrans, tdat_mrg_filt, rms_err] = calc_valerr_hr_vf(pgdat,tdat_mrg,colobj_mrg,1, file_ext);
[tdat_mrg_fullset, rms_fullset] = find_missing_elec_labels_err_vf(nomdat,pgdat, tdat_mrg_filt,colobj_mrg,0, subj_name, file_ext);
%nomps0 = calc_valerr_v2(nomdat,tdat_mrg_fullset,colobj_mrg,1,[scnpth_raw,'/',subj_name]);


% Save the merged data
el = -90*pi/180;% str2double(get(handles.ini_rx,'String'))*pi/180;
Rx        = [1 0 0; 0 cos(el) -sin(el); 0 sin(el) cos(el)];
locs      = double(colobj_mrg.Location);
tdat_mrg  = (Rx*(tdat_mrg'))';
locs      = (Rx*(locs'))';
colobj_mrg= pointCloud(locs, 'Color', colobj_mrg.Color);
tdat      = tdat_mrg;
colobj    = colobj_mrg;
err       = rms_err;
err_fullset= rms_fullset;
newtdats  = [];
file_ini  = ['mergedscan_pg_type', file_ext];

 eval(['save ',savepth,'/file_ini colobj', ...
         ' tdat err err_fullset'])
end
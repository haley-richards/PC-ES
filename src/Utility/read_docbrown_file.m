function mdl = read_docbrown_file(fname)

%-------------------------------------------------------------------------------
% Open the mesh file
fid = fopen(fname);
tline    = fgets(fid);

%-------------------------------------------------------------------------------
stoploop = 0;

while (ischar(tline)) && (stoploop == 0)
    
    %---------------------------------------------------------------------------
    % Get sensor information
    if isempty( strfind(tline,'<sensors>') ) == 0
        sens_ps = zeros(20000,3);
        sens_mn = zeros(20000,1);
        sens_pt = zeros(20000,3);
        sp = 1;
        while isempty( strfind(tline,'</sensors>') ) == 1
            tline    = fgets(fid);
            if isempty( strfind(tline,'<modelPoint>') ) == 0
                % Read in the x,y,z values
                tline = fgets(fid);
                x     = get_brck_val('x',tline);
                tline = fgets(fid);
                y     = get_brck_val('y',tline);
                tline = fgets(fid);
                z     = get_brck_val('z',tline);
                % Get the model number and point type
                tline = fgets(fid);
                im    = get_brck_val('modelNumber',tline);
                tline = fgets(fid);
                ip    = get_brck_val('pointType',tline);
                %-----
                % Record the data
                sens_ps(sp,:) = [x y z];
                sens_mn(sp)   = im;
                sens_pt(sp)   = ip;
                sp = sp + 1;
            end
        end
    end
    %---------------------------------------------------------------------------
    % Get landmark sensor information
    if isempty( strfind(tline,'<landmarkSensors>') ) == 0
        lmrk_ps = zeros(20000,3);
        lmrk_mn = zeros(20000,1);
        lmrk_pt = zeros(20000,3);
        lp = 1;
        while isempty( strfind(tline,'</landmarkSensors>') ) == 1
            tline    = fgets(fid);
            if isempty( strfind(tline,'<modelPoint>') ) == 0
                % Read in the x,y,z values
                tline = fgets(fid);
                x     = get_brck_val('x',tline);
                tline = fgets(fid);
                y     = get_brck_val('y',tline);
                tline = fgets(fid);
                z     = get_brck_val('z',tline);
                % Get the model number and point type
                tline = fgets(fid);
                im    = get_brck_val('modelNumber',tline);
                tline = fgets(fid);
                ip    = get_brck_val('pointType',tline);
                %-----
                % Record the data
                lmrk_ps(lp,:) = [x y z];
                lmrk_mn(lp)   = im;
                lmrk_pt(lp)   = ip;
                lp = lp + 1;
            end
        end
    end
    %---------------------------------------------------------------------------
    % Get landmark sensor information
    if isempty( strfind(tline,'<connections>') ) == 0
        cnxs = zeros(20000,2);
        ic   = 1;
        while isempty( strfind(tline,'</connections>') ) == 1
            tline = fgets(fid);
            cnxs(ic,:) = get_brck_pairval('smc',tline);
            if isnan(cnxs(ic,1)) == 0
            ic = ic + 1;
            end
        end
        
    end
    
    %---------------------------------------------------------------------------
    % get the next line
    tline    = fgets(fid);
    
end
sens_ps = sens_ps(1:(sp-1),:);
sens_mn = sens_mn(1:(sp-1),:);
sens_pt = sens_pt(1:(sp-1),:);
lmrk_ps = lmrk_ps(1:(lp-1),:);
lmrk_mn = lmrk_mn(1:(lp-1),:);
lmrk_pt = lmrk_pt(1:(lp-1),:);
cnxs    = cnxs(1:(ic-1),:);

mdl.sens_ps = [sens_ps; lmrk_ps];
mdl.sens_mn = sens_mn;
mdl.sens_pt = sens_pt;
mdl.lmrk_ps = lmrk_ps;
mdl.lmrk_mn = lmrk_mn;
mdl.lmrk_pt = lmrk_pt;
mdl.cnxs    = cnxs;


%-------------------------------------------------------------------------------
%
%
%-------------------------------------------------------------------------------
function val = get_brck_val(sstr,tline)

i1    = strfind(tline,['<',sstr,'>']);
i2    = strfind(tline,['</',sstr,'>']);
val   = str2double( tline((i1+2+length(sstr)):(i2-1)) );

%-------------------------------------------------------------------------------
function vals = get_brck_pairval(sstr,tline)

i1    = strfind(tline,['<',sstr,'>']);
i2    = strfind(tline,['</',sstr,'>']);
im    = strfind(tline,',');
vals   = [ ...
    str2double( tline((i1+2+length(sstr)):(im-1)) ) ...
    str2double( tline((im+1):(i2-1)) ) ...
    ];




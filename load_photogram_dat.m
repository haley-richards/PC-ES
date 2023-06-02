%-------------------------------------------------------------------------------
% 
% Load the photogrammetry data
% 
%-------------------------------------------------------------------------------
function pgdat = load_photogram_dat(valpth,valfnm)

fid    = fopen([valpth,'/',valfnm]);
fidpos = zeros(3,3);
elpos  = zeros(257,3);
for n = 1:260            
    tline = fgets(fid);
    newStr = split(tline);
    if n <= 3
        fidlab{n}   = newStr{1};
        fidpos(n,:) = [str2double(newStr{2}) str2double(newStr{3}) str2double(newStr{4})];
    else
        elpos(n-3,:) = [str2double(newStr{2}) str2double(newStr{3}) str2double(newStr{4})];
    end
end
fclose(fid)
pgdat.fidlab = fidlab;
pgdat.fidpos = fidpos;
pgdat.elpos  = elpos;
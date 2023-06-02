function ds = signd_pc(p,colobj_mrg)

locs = colobj_mrg.Location;
N    = size(locs,1);
nvcs = colobj_mrg.Normal;
ds   = zeros(size(p,1),1);
size(p,1)
parfor n = 1:size(p,1)
    [tmp,i] = min(sum( (locs - repmat(p(n,:),N,1)).^2,2));
    ds(n)   = sum((p(n,:)-locs(i,:)).*nvcs(i,:));
end
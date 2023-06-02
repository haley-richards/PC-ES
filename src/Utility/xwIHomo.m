function [p] = xwIHomo(homo_p)
%[p] = xwIHomo(homo_p)
%   converts points in homogenous coordinates back to ones.
%
%   Xiaotian (Dennis) Wu

if abs(mean(homo_p(end,:))-1) < 1e-5
    p = homo_p(1:end-1,:);
elseif abs(mean(homo_p(:,end))-1) <1e-5
    p = homo_p(:,1:end-1);
else
    p = homo_p(1:end-1,:);
%     disp('Warning: last row is not all 1s');
end

end


function [disterr_values] = disterr(A,B)
%[disterr_values] = disterr(A,B)
% assume A and B are 3xn if 3x3

sizA = size(A);
sizB = size(B);

if sizA(1) == 3 && sizB(1) == 3
elseif sizA(2) == 3 && sizB(2) == 3
    A = A';
    B = B';
else
    disp('Check Inputs, has to be 3xn or nx3')
    disterr_values = [];
    return
end

disterr_values = sqrt(sum((A-B).^2));

end


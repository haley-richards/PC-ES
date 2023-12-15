function [T, rmserr] = transform_loc(Pa, Pb, flag)
% [T]=transform_loc(Pa,Pb,flag) returns the coordinates transformation relating points
% in matrices Pa and Pb which are the same points in two different frames.
%            | x1 x2 ...|
% Pa or Pb = | y1 y2 ...|
%            | z1 z2 ...|
%            | 1  1  ...|
% T returns the transformation matrix
%
% T is such that Pa = T Pb
%
% The 3rd argument is a flag to select the method of computation.
% Flag=0 for least squares solution.  This may result in very distorted
% space for noisy data),
% Flag=1 for SVD algorithm.
%
% Alex Hartov 06/18/99
%             02/10/09 Cleaned up for ENGG 111 example.
% Xiaotian (Dennis) Wu 12/2016 edited reflection case
% Xiaotian (Dennis) Wu 9/2018 included rmserr calc for registration

switch flag
    case 0 	% LSE result requested
        % Direct inversion method
        T=(Pa*Pb')*inv(Pb*Pb');		% first guess: LSE
        
        
        %% SVD method, most reliable & acurate
    case 1 % Non iterative method based on SVD
        
        % Here we operate with 3x1 vectors and 3x3 rotation matrices
        PA=Pa(1:3,:);		% 3x1 coordinates for the first set of points
        PB=Pb(1:3,:);		% 3x1 coordinates for the 2nd set of points
        Pca=mean(PA,2);		% centroid in a
        Pcb=mean(PB,2);		% centroid in b
        
        %         qa=zeros(3,size(PA,2));
        %         qb=zeros(3,size(PB,2));
        %         H=zeros(3,3);
        %         for i=1:size(PA,2)	% compute vectors from centroids and the H matrix
        %             qa(:,i)=PA(:,i)-Pca;
        %             qb(:,i)=PB(:,i)-Pcb;
        %             H=H+qa(:,i)*qb(:,i)';	% Note: outer product results in a 3x3 matrix
        %
        %         end
        qa = PA - repmat(Pca,1,size(PA,2));
        qb = PB - repmat(Pcb,1,size(PA,2));
        H  = qa*(qb');
        
        %----------------------
        [U,S,V]=svd(H);	% Singular Value Decomposition of H
        
        % XDW 2016-12 edit for reflection correction based on
        % Myronenko, Song 2009 (closed-form solution of rotation matrix)
        C = diag(ones(1,length(U)));
        C(end) = det(U*V');
        X=U*C*V';
        
        % 	X=V*U';			% This is the rotation matrix if all is well
        % 	test=det(X);	% check if we have the rotation or a reflection matrix
        % 	if test == -1	% reflection, needs some adjustment
        % 		V(1:3,3)=-V(1:3,3);
        % 		X=V*U'; 	% NOTE error correction, this step was ommitted before
        % 					% AH 2008-03-19
        % 	end
        
        T=Pca-X*Pcb;		% Compute the offset vector between the two frames of reference
        T=[X T; 0 0 0 1];	% compute the 4x4 affine matrix relating the two frames of references
        
    otherwise
        disp('Use the flag to select between LS (0) and SVD (1) methods');
        return
end

TPb    = T*Pb;
% errs   = disterr(xwIHomo(T*Pb),xwIHomo(Pa));
errs   = disterr(TPb(1:3,:),Pa(1:3,:));
rmserr = sqrt(1/length(errs)*sum(errs.^2));
% rmserr = rms(disterr(xwIHomo(T*Pb),xwIHomo(Pa)));

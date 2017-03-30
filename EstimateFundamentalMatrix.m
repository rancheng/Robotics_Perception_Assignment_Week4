function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2
N = length(x1)
A = zeros(N,9);

for i=1:N
    A(i,:) = [x2(i,1)*x1(i,1) x2(i,1)*x1(i,2) x2(i,1) x2(i,2)*x1(i,1) x2(i,2)*x1(i,2) x2(i,2) x1(i,1) x1(i,2) 1];
end

[~, ~, V] = svd(A);
h = V(:,end);
f = reshape(h,3,3)';
[U,D,V] = svd(f);
D(3,3) = 0;
F = U*D*V';


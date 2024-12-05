
function x = OMP(y, A, k)
% OMP算法
% A: 测量矩阵
% y: 测量结果
% k: 稀疏度
% x: 稀疏信号估计

n = size(A, 2);
r = y;
idx = [];
for i = 1:k
    proj = abs(A'*r);
    [val, pos] = max(proj);
    idx = [idx pos];
    x_hat = pinv(A(:,idx))*y;
    r = y - A(:,idx)*x_hat;
end
x = zeros(n,1);
x(idx) = x_hat;
end
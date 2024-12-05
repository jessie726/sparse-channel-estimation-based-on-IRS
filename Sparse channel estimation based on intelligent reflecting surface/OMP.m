
function x = OMP(y, A, k)
% OMP�㷨
% A: ��������
% y: �������
% k: ϡ���
% x: ϡ���źŹ���

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
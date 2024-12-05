function x_hat = LAOMP(y, A, K, lambda)
% LAOMP算法，带拉格朗日乘子的正交匹配追踪算法

[N, M] = size(A);
r = y;
idx = [];

while length(idx) < K
    % 计算投影
    p = abs(A' * r);
    
    % 选择最大投影对应的列向量
    [~, pos] = max(p);
    
    % 加入索引集合
    idx = [idx, pos];
    
    % 计算线性方程组的解
    if length(idx) > 1
        A_sub = A(:, idx);
        x_hat_sub = pinv(A_sub' * A_sub + lambda * eye(length(idx))) * A_sub' * y;
        x_hat = zeros(M, 1);
        x_hat(idx) = x_hat_sub;
    else
        x_hat = pinv(A(:, idx)) * y;
    end
    
    % 更新残差
    r = y - A * x_hat;
end

end

function x_hat = LAOMP(y, A, K, lambda)
% LAOMP�㷨�����������ճ��ӵ�����ƥ��׷���㷨

[N, M] = size(A);
r = y;
idx = [];

while length(idx) < K
    % ����ͶӰ
    p = abs(A' * r);
    
    % ѡ�����ͶӰ��Ӧ��������
    [~, pos] = max(p);
    
    % ������������
    idx = [idx, pos];
    
    % �������Է�����Ľ�
    if length(idx) > 1
        A_sub = A(:, idx);
        x_hat_sub = pinv(A_sub' * A_sub + lambda * eye(length(idx))) * A_sub' * y;
        x_hat = zeros(M, 1);
        x_hat(idx) = x_hat_sub;
    else
        x_hat = pinv(A(:, idx)) * y;
    end
    
    % ���²в�
    r = y - A * x_hat;
end

end

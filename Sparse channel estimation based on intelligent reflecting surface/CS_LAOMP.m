function [ theta ] = CS_LAOMP( y,A,t,L )
    [M,N] = size(A);
    theta = zeros(N,1);
    At = zeros(M,t);
    Pos_theta = zeros(1,t);
    r_n = y;
    for ii=1:t
        sim = A'*r_n;
        [val1,pos1]=sort(abs(sim),'descend');%pos���������±�
        j = pos1(1:L);
        [jj,ij,io] = intersect(j,Pos_theta);
        j(ij) = [];
        rs = zeros(1,length(j));
        pos_theta = Pos_theta(:,1:ii);
        for i = 1 : length(j)
            dex = j(i);
            rs(i) = norm(LAR_f(y,A,pos_theta,dex,t));
        end
        [val,l] = min(rs);%l�������������д洢fi��ѡ��������
        pos = pos1(l);
        At(:,ii) = A(:,pos);
        Pos_theta(ii) = pos;
        theta_ls = inv(At(:,1:ii)'*At(:,1:ii))*At(:,1:ii)'*y;%��С���˽��α���
        r_n = y - At(:,1:ii)*theta_ls;
    end
    theta(Pos_theta)=theta_ls;
end


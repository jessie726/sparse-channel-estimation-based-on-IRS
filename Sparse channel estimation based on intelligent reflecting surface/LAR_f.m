function [rs] = LAR_f(y,fi,Oinitial,Oupdata,k)
%Oupdata������L���е�һ�������ѭ��ÿ�δ�һ��L�е�һ���������
ss = length(Oinitial);
Oinitial(ss) = Oupdata;%�Ͳ���ʼ�±��������Ϻ���ѡ����L������±�
Onew = Oinitial;
s = length(Onew);%���±꼯�ϳ���
[M,N] = size(fi);
fio = zeros(M,k);
sss = length(Onew);
fio(:,1:sss) = fi(:,Onew);
x = (fio(:,1:sss)'*fio(:,1:sss))^(-1)*fio(:,1:sss)'*y;%pinv(fio)��α��
r_n = y - fio(:,1:sss)*x;%��ʼ�в�
    for i = (s + 1) : k
        product = fi'*r_n;
        [val,pos] = max(abs(product));
        fio(:,i) = fi(:,pos);%�洢��һ��
        Onew(i) = pos;%�洢��һ�е����
        theta_ls = (fio(:,1:i)'*fio(:,1:i))^(-1)*fio(:,1:i)'*y;
        r_n = y - fio(:,1:i)*theta_ls;
    end
    rs = norm(r_n);
end

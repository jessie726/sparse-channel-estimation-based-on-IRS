% --------UPA������������----------
% aziΪ��λ��
% eleΪ����?
% n1��n2�ֱ�Ϊƽ�泤�Ϳ�

function y = steervector(azi,ele,n1,n2)
    N = n1*n2;
    n1 = linspace(0,n1-1,n1);
    n2 = linspace(0,n2-1,n2);
       vec1 = exp(-1j*pi*sin(azi)*cos(ele)*n1);
       vec2 = exp(-1j*pi*sin(ele)*n2);
       y = kron(vec1,vec2);
       y = y.'/sqrt(N);
end


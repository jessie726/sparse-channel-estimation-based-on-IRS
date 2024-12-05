% --------UPA方向向量函数----------
% azi为方位角
% ele为仰角?
% n1和n2分别为平面长和宽

function y = steervector(azi,ele,n1,n2)
    N = n1*n2;
    n1 = linspace(0,n1-1,n1);
    n2 = linspace(0,n2-1,n2);
       vec1 = exp(-1j*pi*sin(azi)*cos(ele)*n1);
       vec2 = exp(-1j*pi*sin(ele)*n2);
       y = kron(vec1,vec2);
       y = y.'/sqrt(N);
end


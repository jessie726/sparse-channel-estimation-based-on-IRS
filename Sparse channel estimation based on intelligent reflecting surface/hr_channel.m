function [hK] = hr_channel(N1,N2,K,Lc,L)


N=N1*N2;
d=0.5;
hK=zeros(N,K);

%%%% generate the common physical angles

N1index=[-(N1-1)/2:1:(N1/2)]'*(2/N1);
N2index=[-(N2-1)/2:1:(N2/2)]'*(2/N2);
index=randperm(N); 
x=ceil(index(1:Lc)/N2);
y=index(1:Lc)-N2*(x-1);
phi1c=N1index(x);
phi2c=N2index(y);

%%%% generate the channel h_{r,k} for each user k

for k=1:K
    alpha = zeros(L,1);
    alpha(1:L) = (normrnd(0, 1, L, 1) + 1i*normrnd(0, 1, L, 1)) / sqrt(2);
    hr=zeros(N,1);

    phi1(1:Lc)=phi1c;
    phi2(1:Lc)=phi2c;
    index=randperm(N);
    x=ceil(index(1:L-Lc)/N2);
    y=index(1:L-Lc)-N2*(x-1);       
    phi1(Lc+1:L)=N1index(x);
    phi2(Lc+1:L)=N2index(y);
    
    for l = 1:L
        a1 = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi1(l));
        a2 = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi2(l));
        a=kron(a1,a2);
        hr = hr + alpha(l)*a;
    end
    
    hK(:,k)= sqrt(N/L)*hr; 
end


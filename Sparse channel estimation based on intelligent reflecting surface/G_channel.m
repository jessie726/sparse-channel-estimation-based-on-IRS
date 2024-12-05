function [G] = G_channel(M1,M2,N1,N2,L)

M=M1*M2;
N=N1*N2;

d=0.5;
G=zeros(M,N);

%%%% generate the physicial angles

M1index=[-(M1-1)/2:1:(M1/2)]'*(2/M1); %原代码把M1错写为N1
M2index=[-(M2-1)/2:1:(M2/2)]'*(2/M2); %原代码把M2错写为N2
N1index=[-(N1-1)/2:1:(N1/2)]'*(2/N1);
N2index=[-(N2-1)/2:1:(N2/2)]'*(2/N2);

index=randperm(N); 
x=ceil(index(1:L)/N2);
y=index(1:L)-N2*(x-1);
phi1=N1index(x);
phi2=N2index(y);

index=randperm(M); 
x=ceil(index(1:L)/M2);
y=index(1:L)-M2*(x-1);
psi1=M1index(x);
psi2=M2index(y);


%%%% generate the gains
alpha = zeros(L,1);
alpha(1:L) = (normrnd(0, 1, L, 1) + 1i*normrnd(0, 1, L, 1)) / sqrt(2);

%%%% generate G channel M*N
for l = 1:L
    
    a1 = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi1(l));
    a2 = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi2(l));
    a=kron(a1,a2);
    b1 = 1/sqrt(M1)*exp(-1i*2*pi*[0:M1-1]'*d*psi1(l));
    b2 = 1/sqrt(M2)*exp(-1i*2*pi*[0:M2-1]'*d*psi2(l));
    b=kron(b1,b2);

    G = G + alpha(l)*b*a.';
end 

G=sqrt(M*N/L)*G;

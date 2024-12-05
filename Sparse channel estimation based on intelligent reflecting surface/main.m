clear all;
close all;
clc;
rng(15000);
SNR_dB = 0:4:20;            %SNR������
SNR_length = length(SNR_dB);%SNR����
ite = 10;                  %��������

N = 81;                     %�����浥Ԫ����?
N1 = 9;
N2 = 9;                     %���䵥Ԫ������N=N1��N2
K  = 1;                     %�������û���
M =  25;                     %��վ������
M1 = 5; 
M2 = 5;                      %M = M1��M2
Q = 100;
s = 1;
P = 2;  %�źŹ���

L_G = 3;          %RIS��BS֮��·��������
L_r = 4;          %user��RIS֮��·��������
alpha_l1 =  (randn(L_G,1)+1j*randn(L_G,1))/sqrt(2);     %RIS-BS���ŵ�����
alpha_r_l2 = (randn(L_r,1)+1i*randn(L_r,1))/sqrt(2);    %user-RIS���ŵ�����

int_num = 128;                                             %������
azimuth_Gr_l1 = randi([0 int_num],1,L_G)*(2*pi)/int_num;   %BS��λ��
elevation_Gr_l1 = randi([0 int_num],1,L_G)*(2*pi)/int_num; %BS����
azimuth_Gt_l1 = randi([0 int_num],1,L_G)*(2*pi)/int_num;   %RIS���䷽λ��
elevation_Gt_l1 = randi([0 int_num],1,L_G)*(2*pi)/int_num; %RIS�������ǧ�
azimuth_r_l2 = randi([0 int_num],1,L_r)*(2*pi)/int_num;    %RIS���շ�λ��
elevation_r_l2 = randi([0 int_num],1,L_r)*(2*pi)/int_num;  %RIS�������ǧ�
%azimuth_Gr_l1 = pi*rand(L_G,2)-pi/2;   %BS��λ��
%elevation_Gr_l1 = pi*rand(L_G,2)-pi/2;  %BS����
%azimuth_Gt_l1 =  pi*rand(L_G,2)-pi/2;  %RIS���䷽λ��?
%elevation_Gt_l1 =  pi*rand(L_G,2)-pi/2;  %RIS�������ǧ�
%azimuth_r_l2 =  pi*rand(L_r,2)-pi/2;    %RIS���շ�λ��?
%elevation_r_l2 =  pi*rand(L_r,2)-pi/2;  %RIS�������ǧ�
%d_BR = 10;  %BS��RIS�ļ��?
%d_RU = 100; %RIS��user�ļ��?
%lamda = 60*10^9; %�źŲ�����d=lamda/2

G = zeros(M,N);    %RIS-BS�ŵ� M��N
h_r = zeros(N,K);  %user-RIS�ŵ� N��K
phi_Q = zeros(N,Q);%Q������������� N��Q
H_k_sparse = zeros(M,N);%H_k�任���ϡ�����?
H_k_sparse_CS = zeros(M,N);
H_k_sparse_CS1 = zeros(M,N);
H_k_sparse_CS2 = zeros(M,N);
H_k_sparse_CS3 = zeros(M,N);
Y_k = zeros(M,Q);    %�۲����Y
W_k = zeros(M,Q);    %��������W

NMSE0 = zeros(SNR_length,1);
NMSE1 = zeros(SNR_length,1);
NMSE2 = zeros(SNR_length,1);
NMSE3 = zeros(SNR_length,1);


%for t1 = 1:M
 %   for t2 = 1:M
 %   D_M(t1,t2)= 1/sqrt(M)*exp(-1j*2*pi*(t1-1)'*(t2-1)/M); %����M��DFT����
  %  end
%end


%for t3 = 1:N
 %   for t4 = 1:N
  %  D_N(t3,t4)= 1/sqrt(N)*exp(-1j*2*pi*(t3-1)'*(t4-1)/N); %����N��DFT����
  %  end
%end
d=0.5;
UN1=(1/sqrt(N1))*exp(-1i*2*pi*[0:N1-1]'*d*[0:N1-1]*(2/N1));
UN2=(1/sqrt(N2))*exp(-1i*2*pi*[0:N2-1]'*d*[0:N2-1]*(2/N2));
D_N=kron(UN1,UN2);

UM1=(1/sqrt(M1))*exp(-1i*2*pi*[0:M1-1]'*d*[-(M1-1)/2:1:(M1/2)]*(2/M1));
UM2=(1/sqrt(M2))*exp(-1i*2*pi*[0:M2-1]'*d*[-(M2-1)/2:1:(M2/2)]*(2/M2));
D_M=kron(UM1,UM2);


%-------test--�ֵ������������֤---------
% D_1 = D_N'*D_N;
% D_2 = D_M'*D_M;
% figure('name', 'D1����ͼ��');
% imagesc(abs(D_1));
% colorbar;
% hold on
% figure('name', 'D2����ͼ��');
% imagesc(abs(D_2));
% colorbar;
%------------------------------------
   for l_one = 1:L_G     %�����ŵ�G

   b1 = steervector(azimuth_Gr_l1(l_one),elevation_Gr_l1(l_one),M1,M2);
   a1 = steervector(azimuth_Gt_l1(l_one),elevation_Gt_l1(l_one),N1,N2);
   G = G + sqrt(M*N/L_G)*alpha_l1(l_one)*b1*a1.';
    
   end

   for l_two = 1:L_r     %�����ŵ�hr
        
   a2 = steervector(azimuth_r_l2(l_two),elevation_r_l2(l_two),N1,N2);
   h_r = h_r + sqrt(1/N)*alpha_r_l2(l_two)*a2;
  
   end

 % ----------------ԭ����H_kͼ��------------------------
    H_k = G*diag(h_r);
    H_k_abs = abs(H_k);
    figure('name', 'H_k����ͼ��');
    imagesc(H_k_abs);
    colorbar;
    set(gca,'YDir','normal');

 %-----------------ԭ�ŵ�H_k��ϡ�����-----------------
    H_k_sparse = D_M'*H_k*(D_N.')'; 
    H_k_sparse_abs = abs(H_k_sparse);

 % ----------------ϡ�����H_k_sparseͼ��--------------
    figure('name', 'H_k_sparse����ͼ��');
    imagesc(H_k_sparse_abs)
    colorbar;
    set(gca,'YDir','normal');
    %------------------------test1---------------------
    
     %for mm = 1:M
     %    for nn = 1:N
      %       if H_k_sparse_abs(mm,nn) < 0.001 %������ֵ?
       %          H_k_sparse_abs(mm,nn) = 0;
        %     else
         %        H_k_sparse_abs(mm,nn) = 1;
          %   end
         %end
     %end
for SNR_index = 1 : SNR_length  
    fprintf("SNR_index = %d\n",SNR_index); 

for iter = 1 : ite
    for q = 1:Q   %----------------Q��ʱ϶------------------
        noise = (randn(M,1)+1j*randn(M,1))/sqrt(2);        %������˹������
        noisepower  = P *10.^(-SNR_dB(SNR_index)/10);      %��������
        wgn_k = sqrt(noisepower)*noise;                    %����
        %-------------
        %for nn = 1:N
         %   phi_Q(nn,q) = exp(-1i*2*pi*(nn-1)*(q-1)/N); %IRS���䵥Ԫ��λʸ����?
        %end
       phi_Q=((rand(N,Q)>0.5)*2-1)/sqrt(N);
        %------------
        %y_k_pre = G*diag(phi_Q(:,q))*h_r*s + wgn_k; %��ʽ�任ǰ�Ľ����ź�y_K
        y_k_post =G*diag(h_r)*phi_Q(:,q)*s + wgn_k;  %��ʽ�任��Ľ����ź�y_K
        Y_k(:,Q) = y_k_post;
        W_k(:,Q) = wgn_k;
        % ------------test2----------------
         %if y_k_post == y_k_pre
          %   fprintf('equal!');
         %end
        
  end
      
       %----------------�ŵ�ģ��----------------
        Y_k = H_k*phi_Q*s+ W_k;                %�ŵ�ģ�͵����Է���
       %--------------ѹ����֪ģ��---------------
        Y_k_CS = (D_M'*Y_k)';                  %ѹ����֪���������
        W_k_CS = (D_M'*W_k)';                  %ѹ����֪����������
        phi_Q_CS = phi_Q'*D_N;                 %ѹ����֪�ĸ�֪����
        %Y_k_CS = phi_Q_CS*H_k_sparse' + W_k_CS; %ѹ����֪�Ļ���ģ��?  
        %-----------------OMP�㷨------------
        H_k_sparse_CS = H_k_sparse_CS.';
        for mm = 1:M
        H_k_sparse_CS(:,mm) =OMP(Y_k_CS(:,mm),phi_Q_CS,L_G*L_r);
        end
        H_k_sparse_CS = H_k_sparse_CS.';
        est_error = H_k_sparse_CS - H_k_sparse;
 
        NMSE0(SNR_index) = NMSE0(SNR_index) + sum(sum(sum(abs(est_error).^2)));
         %-----------------LAOMP�㷨------------
        H_k_sparse_CS1 = H_k_sparse_CS1.';
        for kk = 1:M
        H_k_sparse_CS1(:,kk) =LAOMP(Y_k_CS(:,kk),phi_Q_CS,L_G*L_r,0.2);
        end
        H_k_sparse_CS1 = H_k_sparse_CS1.';
       
        est_error1 = H_k_sparse_CS1 - H_k_sparse;

        NMSE1(SNR_index) = NMSE1(SNR_index) + sum(sum(sum(abs(est_error1).^2)));
        %-----------------COSAMP�㷨-----------------
        H_k_sparse_CS2 = H_k_sparse_CS2.';
        for ll = 1:M
        H_k_sparse_CS2(:,ll) =COSAMP(Y_k_CS(:,ll),phi_Q_CS,L_G*L_r);
        end
        H_k_sparse_CS2 = H_k_sparse_CS2.';
       
        est_error2 = H_k_sparse_CS2 - H_k_sparse;

        NMSE2(SNR_index) = NMSE2(SNR_index) + sum(sum(sum(abs(est_error2).^2)));
        %-----------------SP�㷨-----------------
        H_k_sparse_CS3 = H_k_sparse_CS3.';
        for pp = 1:M
        H_k_sparse_CS3(:,pp) = SP(Y_k_CS(:,pp),phi_Q_CS,L_G*L_r);
        end
        H_k_sparse_CS3 = H_k_sparse_CS3.';
       
        est_error3 = H_k_sparse_CS3 - H_k_sparse;
        NMSE3(SNR_index) = NMSE3(SNR_index) + sum(sum(sum(abs(est_error3).^2)));
end
end
NMSE0 = 10*log10((NMSE0) / ite /10/sum(sum(sum(abs(H_k_sparse).^2))));
NMSE1 = 10*log10((NMSE1)/ ite /10/sum(sum(sum(abs(H_k_sparse).^2))));
NMSE2 = 10*log10((NMSE2) / ite /10/sum(sum(sum(abs(H_k_sparse).^2))));
NMSE3 = 10*log10((NMSE3) / ite/10 /sum(sum(sum(abs(H_k_sparse).^2))));

save testone.mat
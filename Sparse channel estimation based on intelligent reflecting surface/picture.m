clear all;
clc;
load testone.mat;
load testone1.mat;
%plot(SNR_dB,NMSE0,'-or','linewidth',1);
%hold on;
%plot(SNR_dB,NMSE1,'-pb','linewidth',1);
%hold on;
%plot(SNR_dB,NMSE2,'-sm','linewidth',1);
%hold on;
%plot(SNR_dB,NMSE3,'-^c','linewidth',1);
%h = plot(SNR_dB, NMSE0,'-or',SNR_dB, NMSE1,'-pb',SNR_dB,NMSE2,'-sm',SNR_dB,NMSE3,'-^c','linewidth',1);
h = plot(SNR_dB, NMSE0,'-or',SNR_dB, NMSE1,'-pb',SNR_dB,NMSE2,'-sm',SNR_dB,NMSE3,'-^c','linewidth',1);
hold on;
plot(SNR_dB1, NMSE_0,'-.ok',SNR_dB1, NMSE_1,'-.pk',SNR_dB1,NMSE_2,'-.sk',SNR_dB1,NMSE_3,'-.^k');
xlabel('SNR(dB)','FontSize',12,'FontName','Times New Roman');
ylabel('NMSE(dB)','FontSize',12,'FontName','Times New Roman');
%legend('NMSE0-SNR(dB)曲线(OMP)','NMSE1-SNR(dB)曲线(LAOMP)','NMSE2-SNR(dB)曲线(COSAMP)','NMSE3-SNR(dB)曲线(SP)');
legend('NMSE0-SNR(dB)曲线(OMP),K=12','NMSE1-SNR(dB)曲线(LAOMP),K=12','NMSE2-SNR(dB)曲线(COSAMP),K=12','NMSE3-SNR(dB)曲线(SP),K=12','NMSE0-SNR(dB)曲线(OMP),K=20','NMSE1-SNR(dB)曲线(LAOMP),K=20','NMSE2-SNR(dB)曲线(COSAMP),K=20','NMSE3-SNR(dB)曲线(SP),K=20');
%title('算法误差对比图(N=81,M=25,K=12)')
set(gca,'Xlim',[0,20]);
grid on;
hold off;

 %figure(2);
 %semilogy(SNR_dB, MSE,'-ob');
 %xlabel('SNR(dB) ','FontSize',12);
 %ylabel('MSE2 ','FontSize',12);
 %legend('MSE2-SNR(dB)曲线');
 %set(gca,'Xlim',[0,20]);
 %grid on
 %title('RIS channel estimation');


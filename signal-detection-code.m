clear;
clc;
close all;
load('Data.mat');
s_square = sum (s.*s);
for i=1:1:30
    SNR(i) = 10*log10( (A(i)^2) * s_square / sigmaw2);
end
Pfa = [0.1, 0.01, 0.001, 0.0001, 0.00001];

%Calculate theoretical Pd-SNR curve under 5 Pfa values
for i=1:1:30
    for j=1:1:5
        Pd_1(i,j) = qfunc( qfuncinv( Pfa(j) ) -  A(i)*sqrt( s_square / sigmaw2));
        Pd_3(i,j) = qfunc( qfuncinv( Pfa(j)/2 )-  A(i)*sqrt( s_square / sigmaw2))+qfunc( qfuncinv( Pfa(j)/2 )+ A(i) * sqrt(s_square / sigmaw2));
    end 
end


%Calculate numerical Pd-SNR curve under 5 Pfa values
for i=1:1:5
    lambda_1(i) = qfuncinv(Pfa(i))*sqrt(sigmaw2*s_square);
    lambda_2(i) = qfuncinv(Pfa(i)/2)*sqrt(sigmaw2*s_square);
end
sum_mask = sum(mask);
sum_T2=zeros(5,30,10000);
sum_T4=zeros(5,30,10000);
for i=1:1:30
    for j=1:1:10000
        T=0;
        for l=1:1:10
            T=T+s(l)*x(j,l,i);
        end
        for k=1:1:5
            if T>lambda_1(k)
                sum_T2(k,i,j)=1;
            end
            if abs(T)>lambda_2(k)
                sum_T4(k,i,j)=1;
            end
        end
    end
    for m=1:1:5
        A_sum2=0;
        A_sum4=0;
        for n=1:1:10000
            A_sum2=A_sum2+sum_T2(m,i,n)*mask(n);
            A_sum4=A_sum4+sum_T4(m,i,n)*mask(n);
        end
        Pd_2(i,m)=A_sum2/sum_mask;
        Pd_4(i,m)=A_sum4/sum_mask;
    end
end

figure(1);
plot(SNR,Pd_1(:,1),SNR,Pd_1(:,2),SNR,Pd_1(:,3),SNR,Pd_1(:,4),SNR,Pd_1(:,5),'Linewidth',1.5);
hold on;
plot(SNR,Pd_2(:,1),'k--o',SNR,Pd_2(:,2),'k--o',SNR,Pd_2(:,3),'k--o',SNR,Pd_2(:,4),'k--o',SNR,Pd_2(:,5),'k--o','Linewidth',0.9);
title('Theoretical and numerical detection Performance for known A','FontSize',12)
xlabel('Energy-to-noise ratio(dB)','FontSize',10);
ylabel('Probability of detection, P_D','FontSize',10);
legend('P_{FA}=10^{-1}','P_{FA}=10^{-2}','P_{FA}=10^{-3}','P_{FA}=10^{-4}','P_{FA}=10^{-5}','numerical P_D');
axis([min(SNR) max(SNR) 0 1]);
grid on;


figure(2);
plot(SNR,Pd_3(:,1),SNR,Pd_3(:,2),SNR,Pd_3(:,3),SNR,Pd_3(:,4),SNR,Pd_3(:,5),'Linewidth',1.5);
hold on;
plot(SNR,Pd_4(:,1),'k--o',SNR,Pd_4(:,2),'k--o',SNR,Pd_4(:,3),'k--o',SNR,Pd_4(:,4),'k--o',SNR,Pd_4(:,5),'k--o','Linewidth',0.9);
title('Theoretical and numerical detection Performance for unknown A','FontSize',12)
xlabel('Energy-to-noise ratio(dB)','FontSize',10);
ylabel('Probability of detection, P_D','FontSize',10);
legend('P_{FA}=10^{-1}','P_{FA}=10^{-2}','P_{FA}=10^{-3}','P_{FA}=10^{-4}','P_{FA}=10^{-5}','numerical P_D');
axis([min(SNR) max(SNR) 0 1]);
grid on;


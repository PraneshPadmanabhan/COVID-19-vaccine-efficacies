clc
clear all
close all

Rn = rand(3500,3);

Data = importdata('simulatedIndividualParameters.xlsx');
Data1 = Data.data;

ln_beta = log10(Data1(:,6));
[f_beta,x_beta] = ecdf(ln_beta);
figure
plot(x_beta,f_beta,'--')
fiB = interp1(f_beta, x_beta, Rn(:,1), 'linear','extrap')
hold on
plot(fiB,Rn(:,1),'.k','markersize',0.1)
title('log_1_0 (beta)')

delta = Data1(:,9);
[f_delta,x_delta] = ecdf(delta);
figure
plot(x_delta,f_delta,'--')
fid = interp1(f_delta, x_delta, Rn(:,2), 'linear','extrap')
hold on
plot(fid,Rn(:,2),'.k','markersize',0.1)
title('delta')

ln_p = log10(Data1(:,10));
[f_p,x_p] = ecdf(ln_p);
figure
plot(x_p,f_p,'--')
fip = interp1(f_p, x_p, Rn(:,3), 'linear','extrap')
hold on
plot(fip,Rn(:,3),'.k','markersize',0.1)
title('log_1_0 (p)')

params = [10.^fiB fid 10.^fip];
xxx
save('params.mat','params')
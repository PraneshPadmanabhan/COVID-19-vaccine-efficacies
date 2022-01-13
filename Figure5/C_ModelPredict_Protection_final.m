tic 
clc; 
clear all; 
% close all;

NT50model = 565; % Model NT scaling factor see the manuscript 

Data = importdata('V_varyFixedparms_sample2_v1_Scale_10.mat'); % load results

TN = 10000; deli = 0.2; ThT = 100;

for i = 1:1:TN
Pred(i,:) = [Data{i}.NT50_it max(Data{i}.Vsave(:,2)) Data{i}.Vmax];
end

idp = isnan(Pred(:,1));
Pred(idp,:)=[];

Pred(:,1) = Pred(:,1)/NT50model;

% % binning Scaled NT50
ct = 0;
for i = -2.5:deli:2.5
ct = ct+1;
idx = find(log10(Pred(:,1))>=i & log10(Pred(:,1))<i+deli);
Pred1(ct,:) = [i+(0.5*deli) length(find(Pred(idx,2)>ThT)) length(idx) (1-(length(find(Pred(idx,2)>ThT))/length(idx)))*100];
% Check1{ct} = log10(Pred(idx,1));
% if size(idx,1)>0
% Check2(ct,:) = [log10(min(Pred(idx,1))) log10(max(Pred(idx,1)))];
% end
Pred(idx,:)=[];
end

% % Confidence interval
for ij = 1:1:size(Pred1,1)
[~,pci] = binofit(Pred1(ij,2),Pred1(ij,2)+Pred1(ij,3));
VE_95CI(ij,:) = (1 - 2*pci)./(1 - pci);
% VE_95CI_p(ij,:) = [(1 - (((100-Pred1(ij,4))/100)*exp(1.96*sqrt((1/Pred1(ij,2))+(1/Pred1(ij,3)))))) (1 - (((100-Pred1(ij,4))/100)*exp(-1.96*sqrt((1/Pred1(ij,2))+(1/Pred1(ij,3))))))]
end
curve1 = VE_95CI(:,2)*100;
curve2 = VE_95CI(:,1)*100;

% % Protection curve
DataE = [10.^Pred1(:,1) Pred1(:,4) curve1 curve2];

% % Plot results
figure
semilogx(DataE(:,1),DataE(:,2),'-k','linewidth',2)
hold on
plot(DataE(:,1),curve1(:,1),'--k',DataE(:,1),curve2(:,1),'--k','linewidth',1)
ylim([-5,105])
xlim([0.1,10])
xlabel('Scaled NT_5_0')
ylabel('Protection (%)')
set(gca,'fontsize', 18)
axis square
 
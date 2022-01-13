function Median_effect_equation

clear all; close all; clc

% % Load data
data = xlsread('Data_BD236.xlsx'); 

% % Remove NaN values (corresponding to fa<1% and fa>99%; see Figure 1 legend)
idx = isnan(data(:,3));
doseresponse=[10.^data(:,1) data(:,3)*100];
doseresponse(idx,:)=[];

x = log10(doseresponse(:,1));
y = log10(doseresponse(:,2)./(100-doseresponse(:,2)));

% % % parameter estimation (b = estimates, bint = 95% CI) 
K = [x'; ones(1,length(x))]
[b, bint, R, Rint, stats]= regress(y,K.')

% % % Show estimates in the command window
m = b(1)
IC50 = 10^(-b(2)/m)

% % % Plot
xpred = -5:0.1:3;
ypred = b(1)*xpred + b(2);

figure
plot(x,y,'.r',xpred,ypred,'-k','MarkerSize',25,'linewidth',2)
xlabel('log_1_0 [BD-236] (\mug ml^-^1)')
ylabel('log_1_0 [f_a/f_u]')
ylim([-3,3])
xlim([-4.3,2.3])
set(gca,'FontSize',18)

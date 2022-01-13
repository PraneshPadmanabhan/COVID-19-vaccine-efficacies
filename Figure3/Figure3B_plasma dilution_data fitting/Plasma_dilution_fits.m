function Plasma_dilution_fits

clear all; close all; clc

% % for loop to fit data of 3 patients
for i0 = 1:1:3

% % Load data
Name = 'Master_plot_plasma.xlsx'; 
data = xlsread(Name,i0); 
[~,sheets] = xlsfinfo(Name);

% % Remove NaN values (corresponding to fa<1% and fa>99%)
idx = isnan(data(:,3));
doseresponse=[10.^data(:,1) data(:,3)*100];
doseresponse(idx,:)=[];

% % beta0 provides intial guesses. Note that IC50 is log transformed 
[~,idx]=min(abs(doseresponse(:,2)-50));
init_guess = doseresponse(idx,1);
beta0 = [1 log10(init_guess)]; 

% % parameter and 95% CI estimation
[betahat,resid,J]=nlinfit(doseresponse(:,1),doseresponse(:,2),@calc,beta0)
betaci = nlparci(betahat,resid,J)

% % Show estimates in command window
NT50_estimate_CI(i0,:) = [betahat(1) betaci(1,1:2)]
n_estimate_CI(i0,:) = [10^betahat(2) 10.^betaci(2,1:2)]

% % Predict and plot
ln_Dp = [0:0.1:5];
Dp = 10.^ln_Dp;
Pred = calc(betahat,Dp);

figure
semilogx(Dp,Pred,'-',doseresponse(:,1),doseresponse(:,2),'.','MarkerSize',35,'linewidth',2);
hold on 

% % the if else loop plots data not considered for fitting
if i0<3
semilogx(10^data(7,1),data(7,2)*100,'or',10^data(8,1),data(8,2)*100,'or','MarkerSize',8,'linewidth',2)
else
semilogx(10^data(1,1),data(1,2)*100,'or','MarkerSize',8,'linewidth',2)
end

text(80000,7,sheets{i0},'FontSize',14)
xlabel('Plasma dilution')
ylabel('Fraction unaffected, f_u (%)')
ylim([0,110])
xlim([10,100000])
set(gca,'xdir','reverse','FontSize',18)
axis square
end

function F=calc(beta,g)
n = beta(1);
NT50 = 10^beta(2);
F = (g.^n ./ (NT50^n + g.^n))*100;
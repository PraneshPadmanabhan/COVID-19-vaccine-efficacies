clc; clear all; close all; 

% % N is the number of NAbs (here it is 10). D0 is the overall NAb concentration
N = 10; 
D0 = 10; 

% % m and IC50 values of 10 NAbs
IC50 = 10.^[0.76; -0.41; -1.15; 0.29; 0.89; -2.66; 0.49; 0.03; -3.1; 1.3]
m = [1.35; 0.47; 0.4; 0.56; 0.31; 1.5; 0.31; 0.39; 1.27; 0.97]
% % Compute concentration of each NAb
Di(1,1:1:N) = D0/N; 

ct = 0;

% % This for loop computes fraction unaffected at each dilution
for ij0 = 0:0.01:5 
    
Dilution = 10^ij0; % % dilution value
ct = ct + 1;

y = 1;
for i=1:1:N
y = y*(1 - (1 / (1 + ((IC50(i)/(Di(i)/Dilution))^m(i)))));
end

Dil(ct,1) = Dilution; 
fu_Bliss(ct,1) = y*100; % % fu_Bliss is the fraction unaffected calculated using Bliss Independence
end

figure
semilogx(Dil,fu_Bliss,'-r')
ylim([0,100])
xlim([1,100000])
xlabel('Plasma dilution')
ylabel('Fraction unaffected, f_u (%)')
set(gca,'xdir','reverse','FontSize',18)


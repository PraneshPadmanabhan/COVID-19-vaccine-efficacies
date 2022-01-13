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

% % Loewe additivity expression
FF=@(epi) [1 - (Di(1)/(Dilution*IC50(1)))*(((1/epi)-1)^(1/m(1))) ...
             - (Di(2)/(Dilution*IC50(2)))*(((1/epi)-1)^(1/m(2)))...
             - (Di(3)/(Dilution*IC50(3)))*(((1/epi)-1)^(1/m(3)))...
             - (Di(4)/(Dilution*IC50(4)))*(((1/epi)-1)^(1/m(4)))...
             - (Di(5)/(Dilution*IC50(5)))*(((1/epi)-1)^(1/m(5)))...
             - (Di(6)/(Dilution*IC50(6)))*(((1/epi)-1)^(1/m(6)))...
             - (Di(7)/(Dilution*IC50(7)))*(((1/epi)-1)^(1/m(7)))...
             - (Di(8)/(Dilution*IC50(8)))*(((1/epi)-1)^(1/m(8)))...
             - (Di(9)/(Dilution*IC50(9)))*(((1/epi)-1)^(1/m(9)))...
             - (Di(10)/(Dilution*IC50(10)))*(((1/epi)-1)^(1/m(10)))];

options = optimoptions('lsqnonlin','Display','none');
[x] = lsqnonlin(FF,[0.1],[0],[1],options); % % x is fraction affected

Dil(ct,1) = Dilution; 
fu_Loewe(ct,1) = (1-x)*100; % % fu_Loewe is fraction unaffected calculated using Loewe

end

figure
semilogx(Dil,fu_Loewe,'-r')
ylim([0,100])
xlim([1,100000])
xlabel('Plasma dilution')
ylabel('Fraction unaffected, f_u (%)')
set(gca,'xdir','reverse','FontSize',18)


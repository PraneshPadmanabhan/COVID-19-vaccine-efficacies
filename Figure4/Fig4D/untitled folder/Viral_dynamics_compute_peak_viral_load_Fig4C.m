function Viral_dynamics_compute_peak_viral_load_Fig4C
clear; clc; close all
warning('off')

% % Calculations for D0 = 0.1 \mug / ml alone is shown

% % Calculation efficacies of 500 virutal patients-----------------------------

D0 = 0.01; % % D0 is the overall NAb concentration
N = 10; % % N is the number of NAbs
Data = importdata('m_IC50_2.mat'); % % m and IC50 values sampled from the landscape
IC50_save = Data.IC50_save; m_save = Data.m_save;
ct = 0;
IC50s = []; ms = [];

% % 500 represents the number of virtual patient
for r1 = 1:1:500
IC50 = IC50_save(:,r1)'; m = m_save(:,r1)'; Di(1,1:1:N) = D0/N; Dilution = 1; 

IC50s = [IC50s; IC50]; ms = [ms; m];

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
[x] = lsqnonlin(FF,[0.1],[0],[1],options);

ct = ct+1;
Dil(ct,1) = Dilution; 
epi_CI(ct,1) = x;

clear IC50 m
end

epi_CI = [0; epi_CI];
% % epi_CI of 500 virtual patients. 

% -------------------------------------------------------------------------
tspan = [0:0.01:40];

Datasave = [tspan'];

% figure
for i = 1:1:length(epi_CI)
    
% Parameters (rand is used to vary the viral dynamics parameters)----------
p.B = 10^-(7.22185+(0.8*(rand-0.5))); 
p.KX = 4+(2*(rand-0.5));
p.del = 0.6+(1*(rand-0.5));
p.p = 390+(200*(rand-0.5)); 
p.c = 20+(10*(rand-0.5));
p.sX = 1+(1*(rand-0.5));
p.fiX = 100+(200*(rand-0.5));
p.dX = 0.2+(0.1*(rand-0.5));    
p.T0 = 3*10^7;
p.I0 = 1;
p.V0 = ((p.p/p.c)*p.I0);
p.X0 = 0;
p.N = 10;
    
p.epi = epi_CI(i);
% N0 = [p.T0 1 (p.p/p.c) 10^-6 0];
N0 = [p.T0 p.I0 p.V0 p.X0];

[T,N] =  ode45(@(t,y) SARSCoV2(t,y,p),tspan,N0);

Datasave = [Datasave N(:,3)];
psave{i} = p;
assignin('base','V',Datasave);
assignin('base','epi_CI',epi_CI);

% hold on
% % plot(T,(N(:,4)))
% % xlim([0,20]); ylim([0,1])
% plot(T,log10(N(:,3)))
% xlim([0,40]); ylim([0,10])
% xlabel('Time (d)'); ylabel('log_1_0(V)')
% set(gca,'fontsize', 18)

end

save(strcat('V_',num2str(D0),'.mat'),'Datasave','IC50s','ms','psave')

Vpeak = max(Datasave(:,3:end))';
assignin('base','Vpeak',Vpeak); % % Vpeak values used in Fig4C 


function dy = SARSCoV2(~,y,p)

dy = zeros(4,1);

T = y(1);
I = y(2);
V = y(3);
X = y(4);

dy(1) = -p.B*(1-p.epi)*T*V - p.KX*X*T;
dy(2) = p.B*(1-p.epi)*T*V - p.del*I;
dy(3) = p.p*I - p.c*V;
dy(4) = p.sX*(I/(p.fiX+I))*(1-X) - p.dX*X;





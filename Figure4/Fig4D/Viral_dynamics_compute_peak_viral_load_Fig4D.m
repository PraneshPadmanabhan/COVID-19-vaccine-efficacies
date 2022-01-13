function Viral_dynamics_compute_peak_viral_load_Fig4C
clear; clc; close all
warning('off')

% % Calculations for D0 = 0.01 \mug / ml alone is shown

% % Calculation efficacies of 500 virutal patients-----------------------------
Scale = 10; % % The parameter omega (ratio of IC50s in vivo and in vitro)
D0 = 50; % % D0 is the overall NAb concentration
N = 10; % % N is the number of NAbs
Data = importdata('m_IC50_2.mat'); % % m and IC50 values sampled from the landscape
IC50_save = Data.IC50_save; m_save = Data.m_save;
ct = 0;
IC50s_shifted = []; ms = [];

st1=3000; st2=3500; % % st1 and st2 needs to be changed for different D0. 
% % For eg, st1=501; st2=1000; for D0 = 0.1 \mug / ml and so on...

Params_v1 = importdata('params.mat');

% % 500 represents the number of virtual patient
Params_v2 = Params_v1(st1:st2,:);
for r1 = st1:1:st2
IC50 = Scale*IC50_save(:,r1)'; m = m_save(:,r1)'; Di(1,1:1:N) = D0/N; Dilution = 1; 

IC50s_shifted = [IC50s_shifted; IC50]; ms = [ms; m];

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

epi_CI = [epi_CI];
% % epi_CI of 500 virtual patients. 

% -------------------------------------------------------------------------
tspan = [0:0.01:40];

Datasave = [tspan'];

figure
for i = 1:1:length(epi_CI)
    
% Parameters (rand is used to vary the viral dynamics parameters)----------
p.B = Params_v2(i,1); 
p.del = Params_v2(i,2);
p.p = Params_v2(i,3);

p.rhoX = 2+(2*(rand-0.5));
p.k = 4+(4*(rand-0.5));
p.c = 5+(15*(rand));
p.sX = 1+(1*(rand-0.5));
p.fiX = 2+(2*(rand-0.5));
p.dX = 0.2+(0.2*(rand-0.5));

p.T0 = 2666670;
p.I0 = 1/30;
p.V0 = ((p.p/p.c)*p.I0);
p.X0 = 0;
p.N = 10;
    
p.epi = epi_CI(i);
N0 = [p.T0 0 p.I0 p.V0 p.X0];;

[T,N] =  ode23s(@(t,y) SARSCoV2(t,y,p),tspan,N0);
% % or ode45

Datasave = [Datasave N(:,4)];
psave{i} = p;
assignin('base','V',Datasave);
assignin('base','epi_CI',epi_CI);

hold on
% plot(T,(N(:,4)))
% xlim([0,20]); ylim([0,1])
plot(T,log10(N(:,4)))
xlim([0,40]); ylim([0,10])
xlabel('Time (d)'); ylabel('log_1_0(V)')
set(gca,'fontsize', 18)

end

save(strcat('V_',num2str(D0),'.mat'),'Datasave','IC50s_shifted','ms','psave')

Vpeak = max(Datasave(:,2:end))';
assignin('base','Vpeak',Vpeak); % % Vpeak values used in Fig4C 


function dy = SARSCoV2(~,y,p)

dy = zeros(5,1);

T = y(1);
I1 = y(2);
I2 = y(3);
V = y(4);
X = y(5);

dy(1) = -p.B*(1-p.epi)*T*V - p.rhoX*X*T;
dy(2) = p.B*(1-p.epi)*T*V - p.k*I1;
dy(3) = p.k*I1 - p.del*I2;
dy(4) = p.p*I2 - p.c*V;
dy(5) = p.sX*(I2/(I2 + p.fiX))*(1-X) - p.dX*X;


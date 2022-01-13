function B1_scaling_Model_Protection_final
clear; clc; close all; tic

% % This code scales IC50 values from in vitro to in vivo using omega (denoted here as Scale)

tspan = [0:0.01:50];
Totalrun = 10000;
Name = 'V_varyFixedparms_sample2_v1'
DataP = importdata(strcat(Name,'.mat'));
Scale = 10; % % Omega in the paper

for i = 1:1:Totalrun
i
p.NT50_it = DataP{1,i}.NT50_it;
% p.epi = DataP{1,i}.epi;
p.scale = Scale;
p.IC50_it = DataP{1,i}.IC50_it*p.scale;
p.m_it = DataP{1,i}.m_it;
p.A0 = DataP{1,i}.A0;
p.N = DataP{1,i}.N;

IC50 = p.IC50_it; 
m = p.m_it; 
Ai(1,1:1:p.N) = p.A0/p.N; Dilution = 1; 

FF=@(epi) [1 - (Ai(1)/(Dilution*IC50(1)))*(((1/epi)-1)^(1/m(1))) ...
             - (Ai(2)/(Dilution*IC50(2)))*(((1/epi)-1)^(1/m(2)))...
             - (Ai(3)/(Dilution*IC50(3)))*(((1/epi)-1)^(1/m(3)))...
             - (Ai(4)/(Dilution*IC50(4)))*(((1/epi)-1)^(1/m(4)))...
             - (Ai(5)/(Dilution*IC50(5)))*(((1/epi)-1)^(1/m(5)))...
             - (Ai(6)/(Dilution*IC50(6)))*(((1/epi)-1)^(1/m(6)))...
             - (Ai(7)/(Dilution*IC50(7)))*(((1/epi)-1)^(1/m(7)))...
             - (Ai(8)/(Dilution*IC50(8)))*(((1/epi)-1)^(1/m(8)))...
             - (Ai(9)/(Dilution*IC50(9)))*(((1/epi)-1)^(1/m(9)))...
             - (Ai(10)/(Dilution*IC50(10)))*(((1/epi)-1)^(1/m(10)))];

options = optimoptions('lsqnonlin','Display','none');
[x] = lsqnonlin(FF,[0.1],[0],[1],options);
p.epi = x; % Loewe additivity

% Parameters---------------------------------------------------------------
p.B = DataP{1,i}.B;
p.del = DataP{1,i}.del;
p.p = DataP{1,i}.p;

p.T0 = DataP{1,i}.T0;
p.I0 = DataP{1,i}.I0;

p.rhoX = DataP{1,i}.rhoX;
p.k = DataP{1,i}.k;
p.c = DataP{1,i}.c;
p.sX = DataP{1,i}.sX;
p.fiX = DataP{1,i}.fiX;
p.dX = DataP{1,i}.dX;

p.V0 = ((p.p/p.c)*p.I0);
p.M0 = 0;
N0 = [p.T0 0 p.I0 p.V0 p.M0];

[~,N] =  ode45(@(t,y) SARSCoV2(t,y,p),tspan,N0);

p.Vsave = [tspan' N(:,4)];
p.Vmax = max(N(:,4));

ctt = i;
Datasave{ctt} = p;

clearvars -except Datasave tspan Totalrun i DataP VDparam Scale Name
end

save(strcat(Name,'_Scale_',num2str(Scale),'.mat'),'Datasave')

toc
end

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

end

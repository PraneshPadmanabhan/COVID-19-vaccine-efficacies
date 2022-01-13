function B_Model_Protection_final
clear; clc; close all; tic

Rn = rand(10000,3);

Data = importdata('simulatedIndividualParameters.xlsx');
Data1 = Data.data;

ln_beta = log10(Data1(:,6));
[f_beta,x_beta] = ecdf(ln_beta);
figure
plot(x_beta,f_beta,'--')
fiB = interp1(f_beta, x_beta, Rn(:,1), 'linear','extrap')
hold on
plot(fiB,Rn(:,1),'.k','markersize',0.1)
title('beta')

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
title('p')

VDparam = [10.^fiB fid 10.^fip]; % % lines 4-36 generates VD parameters for each individuals

tspan = [0:0.01:50];
Totalrun = 10000;

DataAB = importdata('Ab_new_sample2.mat'); % % NAb landscape sampling

for i = 1:1:Totalrun
i
p.NT50_it = DataAB{1,i};
p.epi = DataAB{2,i}; % Estimated using Loewe additivity
p.IC50_it = DataAB{3,i};
p.m_it = DataAB{4,i};
p.A0 = DataAB{5,i};
p.N = 10;

% Parameters---------------------------------------------------------------
p.B = (VDparam(i,1));
p.del = (VDparam(i,2));
p.p = (VDparam(i,3));

p.T0 = 2666670;
p.I0 = 1/30;

p.rhoX = 2+(2*(rand-0.5));
p.k = 4+(4*(rand-0.5));
% p.c = 10;
p.c = 5+(15*rand);
p.sX = 1+(1*(rand-0.5));
p.fiX = 2+(2*(rand-0.5));
p.dX = 0.2+(0.2*(rand-0.5));

% p.rhoX = 2;
% p.k = 4;
% p.c = 10;
% p.sX = 1;
% p.fiX = 2;
% p.dX = 0.2;

p.V0 = ((p.p/p.c)*p.I0);
p.X0 = 0;
N0 = [p.T0 0 p.I0 p.V0 p.X0];

% [~,N] =  ode23s(@(t,y) SARSCoV2(t,y,p),tspan,N0);
[~,N] =  ode45(@(t,y) SARSCoV2(t,y,p),tspan,N0);

p.Vsave = [tspan' N(:,4)];
p.Vmax = max(N(:,4));

ctt = i;
Datasave{ctt} = p;

clearvars -except Datasave tspan Totalrun i DataAB VDparam
end

save(strcat('V_varyFixedparms_sample2_v1.mat'),'Datasave')

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

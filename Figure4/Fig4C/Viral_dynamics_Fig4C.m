function Viral_dynamics_Fig4C
clear; clc; 
close all

% Efficacies explored------------------------------------------------------ 
epi_CI = [0 0.4 0.6 0.8 0.95 0.98];

% Parameters--------------------------------------------------------------- 

p.B = 2.22*10^-8; 
p.rhoX = 2;
p.del = 0.96;
p.p = 7276.4; 
p.c = 10;
p.k = 4;
p.sX = 1;
p.fiX = 2;
p.dX = 0.2;

% Initial conditions-------------------------------------------------------
p.T0 = 2666670;
p.I0 = 1/30;
p.V0 = (p.p/p.c)*p.I0;
p.M0 = 0;

% -------------------------------------------------------------------------
tspan = [0:0.01:55];
Datasave = [tspan'];

% % prediction viral dynamics for different efficacies
for i = 1:1:length(epi_CI)
p.epi = epi_CI(i);
N0 = [p.T0 0 p.I0 p.V0 p.M0];

[T,N] =  ode23s(@(t,y) SARSCoV2(t,y,p),tspan,N0);

Datasave = [Datasave N(:,4)];
assignin('base','V',Datasave);
assignin('base','epi_CI',epi_CI);
% assignin('base','name',N);

hold on
plot(T,log10(N(:,4)))
xlim([0,35]); ylim([0,8])
xlabel('Time (d)'); ylabel('log_1_0 [V]')
set(gca,'fontsize', 18)
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
dy(5) = p.sX*(I2/(p.fiX+I2))*(1-X) - p.dX*X;







close all; clear all; clc
load seismic.mat
addpath Main_Routines/

%% generate a mesh

% define the polygon we want to mesh
x_min = 0; x_max = 3; 
z_min = 0; z_max = 3; 
pgon = polyshape([x_min x_min x_max x_max], [z_min z_max z_max z_min]);
% pgon = polyshape([0 0 0:0.1:2*pi 2*pi 2*pi], [0 1 0.5*sin(0:0.1:2*pi)+2 1 0]);

% generate a linear mesh
%Hmax = 0.025;
Hmax = 0.04;
%Hmax = 0.08;
Hmin = Hmax/2;
[Elements, Nodes, model, mesh] = GenerateLinearMesh(pgon,Hmin,Hmax);

N_el = length(Elements);
N = length(Nodes);
n = size(Elements,1);  

% plot the mesh
figure(1)
pdemesh(model)

% dfine the velocity model and source locatoin
%c = @(x,z) 1 + x - x + z - z; 
% c = @(x,z) (2.5 - z/2).^2;
c = @(x,z) 1 + x - x ;
figure(100)
pdeplot(model,'XYData',c(Nodes(1,:),Nodes(2,:)),'Mesh','off','ColorMap',jet)

% create f(x,z) forcing function
s = 0.001;
x_s = 1.5;
z_s = 1.5;
f = @(x,z) exp(-((x-x_s).^2)/(2*s^2)).*exp(-((z-z_s).^2)/(2*s^2));

% dfine time and source function 
T_max = 3.5; 
dt = sqrt(1/2)*Hmin/max(c(Nodes(1,:),Nodes(2,:))); 
NT = round(T_max/dt);
t = dt*(0:NT-1);
f_0 = 2.5;
t_0 = 0 + 1/f_0;
f_t = -(1 - 2*pi^2*f_0^2*(t-t_0).^2).*exp(-pi^2*f_0^2*(t-t_0).^2);

% plt the source time function f(t)
figure(10)
plot(t,f_t,'k','LineWidth',2)
xlabel('time')
ylabel('amplitude')
title('Source Time Function')
grid on

%% Assembly

[F, M, K] = Assembly(Elements,Nodes,f,c); 

% Set boundary conditions
% say which edges you want to set the aborbing boundary conditions at 
edges_ABC = [1 2 3 4];
[M_tilda,M,K,F] = Set_BCs(model,mesh,Elements,Nodes,edges_ABC,M,K,F,c);

%%

freq = 10;
L = (-(freq*2*pi)^2*M + 1i*freq*2*pi*M_tilda + K);
tic
F = 0*F; F(round(end/2)) = 1;
C1 = L\F;
F = 0*F; F(round(end)) = 1i;
C2 = L\F;
C = real(C1.*C2);
% C = M\-K*C + C;
toc

%%
figure; 
subplot(3,1,1)
pdeplot(model,'XYData',real(C1),'Mesh','off','ColorMap',jet)
caxis(max(abs(C1))*[-1 1]/4)
subplot(3,1,2)
pdeplot(model,'XYData',real(C2),'Mesh','off','ColorMap',jet)
caxis(max(abs(C2))*[-1 1]/4)
subplot(3,1,3)
pdeplot(model,'XYData',C,'Mesh','off','ColorMap',cmap)
caxis(max(abs(C))*[-1 1]/4)

    
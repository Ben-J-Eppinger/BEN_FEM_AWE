
close all; clear all; clc
load seismic.mat
addpath Main_Routines/

x = 0:0.01:(3/2)*pi;
r = 0.5*(-sin(x) + 1);
b = 0.5*(cos(x) + 1);
g = 0.5*(r+b);%abs(sinc(x-(3/4)*pi));
cmap = [r' g' b'];

%% generate a mesh

% define the polygon we want to mesh
x_min = 0; x_max = 3; 
z_min = 0; z_max = 3; 
pgon = polyshape([x_min x_min x_max x_max], [z_min z_max z_max z_min]);
% pgon = polyshape([0 0 0:0.1:2*pi 2*pi 2*pi], [0 1 sin(0:0.1:2*pi)+2 1 0]);

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
% c = @(x,z) 1 + x - x + z - z; 
% c = @(x,z) (2.5 - z/2).^2;
c = @(x,z) 1 + 0*x + 0*z;

figure(100)
pdeplot(model,'XYData',c(Nodes(1,:),Nodes(2,:)),'Mesh','off','ColorMap',jet)

% create f(x,z) forcing function
s = 0.001;
x_s = 0.1;
z_s = 0.1;
f = @(x,z) exp(-((x-x_s).^2)/(2*s^2)).*exp(-((z-z_s).^2)/(2*s^2));

% dfine time and source function 
T_max = 3.5; 
dt = sqrt(1/2)*Hmin/max(c(Nodes(1,:),Nodes(2,:))); 
NT = round(T_max/dt);
t = dt*(0:NT-1);
f_0 = 2;
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


%% integgrate over time

% parameter to assert which timespetting method to use
Damp = 1; 

% initialize measure of energy over time
E = nan(NT,1);

% initialize degrees of freedom
C = zeros(size(F));
C_old = zeros(size(C));

tic 
count = 1;
p_done = 0;
fprintf('Modeling Wavefield: %6.0f%%',p_done)
for i = 1:NT
    
    if Damp == 0
        C_new = (dt^2)*(M\(f_t(i)*F - K*C)) + 2*C - C_old;
    else
        C_new = (M + dt*M_tilda)\...
            ((dt^2)*(f_t(i)*F - K*C) + M*(2*C - C_old) + dt*M_tilda*C);
    end
    
    C_old = C;
    C = C_new; 
    
    E(i) = sum(C.^2);
    
    p_done = 100*i/NT;
    fprintf('\b\b\b\b\b\b\b\b%3.0f%% || ',p_done)

    if mod(i,floor(NT/9)) == -1%0
        
        if count == 1
           scale = max(abs(C)); 
        end
        
%         p_done = 100*i/NT;
%         fprintf('\b\b\b\b\b\b\b\b%3.0f%% || ',p_done)
        
        % plot a snapshot of the wavefield
        figure(30);
        subplot(3,3,count); count = count + 1; 
        pdeplot(model,'XYData',C,'Mesh','off','ColorMap',seismic)
        %pdeplot(model,'XYData',C,'Mesh','off','ColorMap',cmap)
        caxis([-1 1]*scale*0.5)
        %caxis([-1 1]*10^-8)
        set(gcf, 'Position',  [200, 200, 1200, 400])
        
        % plot E(t)
        figure(40)
        plot(t,E,'k','LineWidth',2)
        grid on
        xlabel('time')
        ylabel('\int_\Omega u^2 dx dz')
        set(gcf, 'Position',  [1200, 400, 1400, 600])
                 
    end
end
toc

    
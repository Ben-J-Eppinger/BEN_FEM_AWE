
close all; clear all; clc
load seismic.mat

%% generate a mesh

% define the polygon we want to mesh
x_min = 0; x_max = 2; 
z_min = 0; z_max = 2; 
pgon = polyshape([x_min x_min x_max x_max], [z_min z_max z_max z_min]);

% generate a linear mesh
Hmin = 0.025;
Hmax = 0.050;
[Elements, Nodes, model, mesh] = GenerateLinearMesh(pgon,Hmin,Hmax);

N_el = length(Elements);
N = length(Nodes);
n = size(Elements,1); 
elements_per_node = Calc_Elements_Per_Node(Elements, Nodes); 

% plot the mesh
figure(1)
pdemesh(model)

%%
vp = @(x,z) 1; 
s = 0.001;
cen = 1;
f = @(x,z) exp(-((x-cen).^2)/(2*s^2)).*exp(-((z-cen).^2)/(2*s^2));

%% Assembly

[F, M, K] = Assembly(Elements,Nodes,model,mesh,f,vp); 

%% plot where the source function is

figure(2);
pdeplot(model,'XYData',F,'Mesh','on','ColorMap',seismic)
caxis([-1 1]*max(abs(F))*1)

%% unitialize time variables and integgrate over time

dt = sqrt(1/2)*Hmin; 
NT = 90;
t = dt*(0:NT);
f_0 = 3;
t_0 = 0 + 1/f_0;
f_t = -(1 - 2*pi^2*f_0^2*(t-t_0).^2).*exp(-pi^2*f_0^2*(t-t_0).^2);

C = zeros(size(F));
[mem_x, mem_z, dx_memx] = deal(zeros(size(C))); 
C_old = zeros(size(C));

tic 
for i = 1:NT
    
    C_new = (dt^2)*(M\(f_t(i)*F - K*C)) + 2*C - C_old;
    C_old = C;
    C = C_new; 

    if mod(i,round(NT/10)) == 0
        
        p_done = 100*i/NT;
        fprintf('Time integration is %f percent done \n',p_done)
        
        figure(30);
        pdeplot(model,'XYData',C,'Mesh','on','ColorMap',seismic)
        caxis([-1 1]*max(abs(C))*1)
        
    end
end
toc

    
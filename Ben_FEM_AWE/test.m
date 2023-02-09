
close all; clear all; clc
load seismic.mat

%% generate a mesh

% define the polygon we want to mesh
x_min = 0; x_max = 2; 
z_min = 0; z_max = 2; 
%pgon = polyshape([0 0 H/2 H H], [0 H 0.75*H H 0]);
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
%s = 0.25;
s = 0.001;
cen = 1;
f = @(x,z) exp(-((x-cen).^2)/(2*s^2)).*exp(-((z-cen).^2)/(2*s^2));

%% Assembly

[F, M, K] = Assembly(Elements,Nodes,model,mesh,f,vp); 

%% plot where the source function is

figure(2);
pdeplot(model,'XYData',F,'Mesh','on','ColorMap',seismic)
caxis([-1 1]*max(abs(F))*1)

%% integgrate over time

dt = sqrt(1/2)*Hmin; 
NT = 90;
t = dt*(0:NT);
f_0 = 3;
t_0 = 0 + 1/f_0;
f_t = -(1 - 2*pi^2*f_0^2*(t-t_0).^2).*exp(-pi^2*f_0^2*(t-t_0).^2);

%%

[a_x,a_z,b_x,b_z] = Generate_CPMLs(x_max,x_min,z_max,z_min,1,dt,f_0,0.4,Nodes,1);
figure(5);
subplot(2,2,1)
pdeplot(model,'XYData',a_x,'Mesh','on','ColorMap',jet)
title('a_x')
subplot(2,2,2)
pdeplot(model,'XYData',a_z,'Mesh','on','ColorMap',jet)
title('a_z')
subplot(2,2,3)
pdeplot(model,'XYData',b_x,'Mesh','on','ColorMap',jet)
title('b_x')
subplot(2,2,4)
pdeplot(model,'XYData',b_z,'Mesh','on','ColorMap',jet)
title('b_z')


%% integgrate over time

C = zeros(size(F));
[mem_x, mem_z, dx_memx] = deal(zeros(size(C))); 
C_old = zeros(size(C));

tic 
for i = 1:NT
    
    
    %C_new = (dt^2)*(M\(f_t(i)*F-K*C)) + 2*C - C_old;
    C_new = (dt^2)*(M\(f_t(i)*F - K*C)) + 2*C - C_old;
    C_new = C_new + a_x.*C_new + a_z.*C_new;
    C_old = C;
    C = C_new; 
    
    [d_x_u, d_z_u] = Calc_Gradient(Elements, Nodes, elements_per_node, C);
    mem_x = (b_x.*mem_x + a_x.*d_x_u);  
    
    [dx_memx, ~] = Calc_Gradient(Elements, Nodes, elements_per_node, mem_x);
    
    

    if mod(i,round(NT/10)) == 0
        
        p_done = 100*i/NT;
        fprintf('Time integration is %f percent done \n',p_done)
        
        figure(30);
        pdeplot(model,'XYData',C,'Mesh','on','ColorMap',seismic)
        caxis([-1 1]*max(abs(C))*1)
        
        figure(40)
        %subplot(1,2,1)
        pdeplot(model,'XYData',dt*dx_memx,'Mesh','on','ColorMap',seismic)
        caxis([-1 1]*max(abs(C))*1)
        %subplot(1,2,2)
        %pdeplot(model,'XYData',d_z_u,'Mesh','on','ColorMap',seismic)
        
    end
end
toc

%% compute gradient of u 

%C = (sin(Nodes(1,:)).*sin(Nodes(2,:)))'; 
% figure(2);
% pdeplot(model,'XYData',C,'Mesh','on','ColorMap',seismic)
% caxis([-1 1]*max(abs(C))*1)

% figure(3);`
% pdeplot(model,'XYData',d_x_u,'Mesh','on','ColorMap',seismic)
%caxis([-1 1]*max(abs(C))*1)
    
%         figure(4)
%         subplot(1,2,1)
%         pdeplot(model,'XYData',d_x_u,'Mesh','on','ColorMap',seismic)
%         subplot(1,2,2)
%         pdeplot(model,'XYData',d_z_u,'Mesh','on','ColorMap',seismic)
%     
%  NID = findNodes(mesh,'region','Edge',1);
%  EID = findElements(mesh,'attached',NID); 
    
    


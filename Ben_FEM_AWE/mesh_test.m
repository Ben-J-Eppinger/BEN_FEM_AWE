close all; clear all; clc

% define the polygon we want to mesh
x_min = 0; x_max = 3; 
z_min = 0; z_max = 3; 
pgon = polyshape([x_min x_min x_max x_max], [z_min z_max z_max z_min]);
%pgon = polyshape([0 0 0:0.1:2*pi 2*pi 2*pi], [0 1 sin(0:0.1:2*pi)+2 1 0]);

% generate a linear mesh
Hmax = 0.3;
Hmin = Hmax/2;
[Elements, Nodes, model, mesh] = GenerateLinearMesh(pgon,Hmin,Hmax);

N_el = length(Elements);
N = length(Nodes);
n = size(Elements,1);  

% plot the mesh
figure(1)
subplot(2,2,1)
pdemesh(model)
%%
clear all; clc

% define the polygon we want to mesh
x_min = 0; x_max = 3; 
z_min = 0; z_max = 3; 
pgon = polyshape([x_min x_min x_max x_max], [z_min z_max z_max z_min]);
%pgon = polyshape([0 0 0:0.1:2*pi 2*pi 2*pi], [0 1 sin(0:0.1:2*pi)+2 1 0]);

% generate a linear mesh
Hmax = 0.15;
Hmin = Hmax/2;
[Elements, Nodes, model, mesh] = GenerateLinearMesh(pgon,Hmin,Hmax);

N_el = length(Elements);
N = length(Nodes);
n = size(Elements,1);  

% plot the mesh
figure(1)
subplot(2,2,2)
pdemesh(model)


%%
clear all; clc

% define the polygon we want to mesh
x_min = 0; x_max = 3; 
z_min = 0; z_max = 3; 
%pgon = polyshape([x_min x_min x_max x_max], [z_min z_max z_max z_min]);
pgon = polyshape([0 0 0:0.1:2*pi 2*pi 2*pi], [0 1 sin(0:0.1:2*pi)+2 1 0]);

% generate a linear mesh
Hmax = 0.3;
Hmin = Hmax/2;
[Elements, Nodes, model, mesh] = GenerateLinearMesh(pgon,Hmin,Hmax);

N_el = length(Elements);
N = length(Nodes);
n = size(Elements,1);  

% plot the mesh
figure(1)
subplot(2,2,3)
pdemesh(model)
%%
clear all; clc

% define the polygon we want to mesh
x_min = 0; x_max = 3; 
z_min = 0; z_max = 3; 
%pgon = polyshape([x_min x_min x_max x_max], [z_min z_max z_max z_min]);
pgon = polyshape([0 0 0:0.1:2*pi 2*pi 2*pi], [0 1 sin(0:0.1:2*pi)+2 1 0]);

% generate a linear mesh
Hmax = 0.15;
Hmin = Hmax/2;
[Elements, Nodes, model, mesh] = GenerateLinearMesh(pgon,Hmin,Hmax);

N_el = length(Elements);
N = length(Nodes);
n = size(Elements,1);  

% plot the mesh
figure(1)
subplot(2,2,4)
pdemesh(model)

function [Elements,Nodes,model,mesh] = GenerateLinearMesh(pgon,Hmin,Hmax)
%GenerateLinearMesh Generates a linear mesh from an input polygon
%   inputs:
%   pgon: Polygon to mesh
%   Hmin: the minimum allowed element size
%   Hmax: the maximum allowed element size
%
%   outputs:
%   Elements: linear elements
%   Nodes: Linear ndoes
%   Model: PDE model (for plotting)
%   mesh: data structure with mesh properties

tr = triangulation(pgon);

model = createpde;

tnodes = tr.Points';
telements = tr.ConnectivityList';

geometryFromMesh(model,tnodes,telements);

mesh = generateMesh(model,'GeometricOrder','linear','Hmin',Hmin,'Hmax',Hmax);

Elements = mesh.Elements;
Nodes = mesh.Nodes;

end

    
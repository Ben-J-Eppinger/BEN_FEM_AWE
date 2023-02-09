function [elements_per_node] = Calc_Elements_Per_Node(Elements,Nodes)
%Calc_Elements_Per_Node Summary of this function goes here
%   Detailed explanation goes here

N = size(Nodes,2);
elements_per_node = zeros(1,N); 

% loop over all nodes
for i = 1:N
   elements_per_node(i) = sum(sum(Elements == i));  
end

end


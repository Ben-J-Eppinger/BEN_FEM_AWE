
function [a_x,a_z,b_x,b_z] = Generate_CPMLs(x_max,x_min,z_max,z_min,max_vp,dt,f_0,L,Nodes,plt_show)

spac = min([x_max-x_min z_max-z_min])/100;
[x, z] = meshgrid(x_min:spac:x_max, z_min:spac:z_max);
n_x = size(x,2);
n_z = size(z,1);

% paramerters for gaussian tamping
xbl = L;

xi = linspace(0,xbl,round(xbl/spac));

alpha = pi*f_0*(1-xi);

R = 10^-5;
d_0 = log(1/R)*3/(2*xbl);
d = d_0*max_vp*(xi.^2);

b = exp(-(d+alpha)*dt);

a = (b-1).*d./(d+alpha);

%

a_z = zeros(n_x,n_z);
b_z = zeros(n_x,n_z);

xbl = length(a);
zbl = xbl;

for i = 1:n_x
    a_z(i,1:zbl) = flip(a)';
    a_z(i,end-zbl+1:end) = a';
    b_z(i,1:zbl) = flip(b)';
    b_z(i,end-zbl+1:end) = b';
end

a_x = zeros(n_x,n_z);
b_x = zeros(n_x,n_z);

for i = 1:n_z
    a_x(1:xbl,i) = flip(a)';
    a_x(end-xbl+1:end,i) = a';
    b_x(1:xbl,i) = flip(b)';
    b_x(end-xbl+1:end,i) = b';
end

%
a_x = a_x';
a_z = a_z';
b_x = b_x';
b_z = b_z';

a_x = scatteredInterpolant(x(:),z(:),a_x(:),'linear');
a_x = a_x(Nodes(1,:), Nodes(2,:))';

a_z = scatteredInterpolant(x(:),z(:),a_z(:),'linear');
a_z = a_z(Nodes(1,:), Nodes(2,:))';

b_x = scatteredInterpolant(x(:),z(:),b_x(:),'linear');
b_x = b_x(Nodes(1,:), Nodes(2,:))';

b_z = scatteredInterpolant(x(:),z(:),b_z(:),'linear');
b_z = b_z(Nodes(1,:), Nodes(2,:))';

if plt_show == 1
    figure;
    subplot(2,2,1)
    plot(a, 'Linewidth', 2)
    grid on
    title('a')
    subplot(2,2,2)
    plot(b, 'Linewidth', 2)
    grid on
    title('b')
    subplot(2,2,3)
    plot(alpha, 'Linewidth', 2)
    grid on
    title('alpha')
    subplot(2,2,4)
    plot(d, 'Linewidth', 2)
    grid on
    title('d')
end

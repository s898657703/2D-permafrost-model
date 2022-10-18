function [ipx,ipy] = FMESH(Mesh)

%%  Generate mesh in horizontal direction
x_mesh = Mesh.ipx;
m   = size(x_mesh,1);
ipx = x_mesh(1,1);

for j = 1:m
    
    dx   = x_mesh(j,3);
    ipx1 = (x_mesh(j,1) + dx : dx : x_mesh(j,2))';
    ipx  = [ipx ; ipx1];
    
end

%%  Generate mesh in vertical direction
y_mesh = Mesh.ipy;
m = size(y_mesh,1);
ipy = y_mesh(1,1);

for j = 1:m
    
    dy   = y_mesh(j,3);
    ipy1 = (y_mesh(j,1) + dy : dy : y_mesh(j,2))';
    ipy  = [ipy ; ipy1];
    
end

end


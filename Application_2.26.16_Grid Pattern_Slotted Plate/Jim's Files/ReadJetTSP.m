clear
clc

load Run2TMesh.mat;
Tdata = Run2MeshStackT;

load XMesh.mat;
load YMesh.mat;
load ZMesh.mat;

% loop through all images in sequence
tmp = size(Tdata);
total = tmp(3);

for k = 1:total,
    img = Tdata(:,:,k);
    figure(1)
    hold off
    imagesc(Xr(1,:),Yr(:,1),img,[60 75]); % Run 2 range
%    imagesc(Xr(1,:),Yr(:,1),img,[50 60]); % Run 3 range
%    imagesc(Xr(1,:),Yr(:,1),img,[45 54]); % Run 4 range
    axis equal
    xlabel('X (inches)');
    colormap jet
    colorbar
    pause(0.5);
end

clear all
clc
close all




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Arvind K.\ Saibaba 
% North Carolina State Unviersity
% If you use this code in any form, please cite "Tensor-based flow 
% reconstruction from optimally located sensor measurements" 
% by M. Farazmand and A. Saibaba, J. Fluid Mech. (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Read the data
udata = ncread('~/Desktop/tangaroa.nc','u') ;
uctrain = double(udata(1:2:end, 1:2:end, 1:2:end,1:150)); %roughly 75%
uctest = double(udata(1:2:end, 1:2:end, 1:2:end,151:end));


%% Perform mean subtraction
ucm = mean(uctrain, 4);
uctrain = bsxfun(@minus, uctrain, ucm);
clear udata


%% Training: Tensor-DEIM
ulst = cell(1,3);
plst = cell(1,3);
wlst = cell(1,3);
ind = [5, 5, 5];
for j = 1:3
    [u,s,~] = svdsketch(tens2mat(uctrain,j), 1.e-2);
    u = u(:,1:ind(j));
    [~,~,p] = qr(u',0);
    ulst{j} = u;
    plst{j} = p(1:ind(j));
    wlst{j} = u/u(plst{j},:);
end


%% Testing error: Tensor DEIM
err = zeros(size(uctest,4),1);
for k = 1:size(uctest,4)
    ten = uctest(:,:,:,k) - ucm;
    tens = ten(plst{1},plst{2},plst{3});
    
    tenr = tmprod(tens, wlst, 1:3);
    err(k) = norm(ten(:)-tenr(:))/norm(ten(:)+ucm(:));
end

%% Visualize flow field
ten = uctest(:,:,:,25) - ucm;     % Arbitrary test snapshot
tens = ten(plst{1},plst{2},plst{3});
tenr = tmprod(tens, wlst, 1:3);
cmax = max([ten(:) + ucm(:); tenr(:) + ucm(:)]);
cmin = min([ten(:) + ucm(:); tenr(:) + ucm(:)]);

x = linspace(-0.35, 0.65, 150);
y = linspace(-0.3, 0.3, 90);
z = linspace(-0.5, -0.3, 60);
[X,Y,Z] = meshgrid(y, x, z);

xslice = 0.0;
yslice = [0.0,0.4];
zslice = -0.35;
figure,

fsize = 18;
rect = [0,0, 14, 4];
set(gcf, 'Units', 'inches');
set(gcf, 'OuterPosition',rect);
set(gcf, 'Position', rect);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'defaultaxesfontsize', fsize);
set(gcf, 'defaulttextfontsize', fsize);
set(0, 'DefaultAxesFontName','Times New Roman');
set(0, 'defaultTextFontName','Times New Roman');

subplot(1,2,1)
slice(X,Y,Z,uctest(:,:,:,1),xslice,yslice,zslice)
hold on
colorbar
colormap('french')
shading interp
[Px,Py,Pz] = ndgrid(y(plst{2}), x(plst{1}), z(plst{3}));
plot3(Px(:),Py(:),Pz(:), 'k.')
xlabel('y')
ylabel('x')
zlabel('z')
set(gca, 'FontSize', 14);
axis equal
caxis([cmin,cmax])
view([-42.2468749966196 6.82274117589402]);
title('True u-velocity', 'FontSize', 16)
subplot(1,2,2)
slice(X,Y,Z,tenr+ucm,xslice,yslice,zslice)
hold on
colorbar
colormap('french')
shading interp
[Px,Py,Pz] = ndgrid(y(plst{2}), x(plst{1}), z(plst{3}));
plot3(Px(:),Py(:),Pz(:), 'k.')
xlabel('y')
ylabel('x')
zlabel('z')
set(gca, 'FontSize', 14);
axis equal
caxis([cmin,cmax])
view([-42.2468749966196 6.82274117589402]);
re = norm(ten(:)-tenr(:))/norm(ten(:)+ucm(:));
str = strcat(strcat('Tensor-DEIM, RE = ',num2str(100*re, '%.2f')), '%');
title(str, 'FontSize', 16)
print -depsc tangaroa





%% This is not included in the paper but maybe useful


%% Visualize: POD-DEIM reconstruction
r = 125;    % Corresponds to 5 x 5 x 5 sensors
S = tens2mat(uctrain,4)';
[u,p] = poddeimrand(S, r);
vec = ten(:);
podr = u*(u(p,:)\vec(p));   % reconstruction
rer = norm(podr(:)-ten(:))/norm(ten(:) + ucm(:));

figure,
slice(X,Y,Z, reshape(podr + ucm(:), size(ten)),xslice,yslice,zslice)
hold on
colorbar
colormap('french')
shading interp
pts = [X(:),Y(:),Z(:)];
plot3(pts(p,2), pts(p,1), pts(p,3), 'k.')
xlabel('y')
ylabel('x')
zlabel('z')
set(gca, 'FontSize', 14);
axis equal
caxis([cmin,cmax])
view([-42.2468749966196 6.82274117589402]);
str = strcat(strcat('POD-DEIM, RE = ',num2str(100*rer, '%.2f')), '%');
title(str, 'FontSize', 16)



%% Testing error: POD-DEIM
for k = 1:size(uctest,4)
    ten = uctest(:,:,:,k)-ucm; vec = ten(:);
    podr = u*(u(p,:)\vec(p));
    rer(k) = norm(ten(:)-podr(:))/norm(ten(:) + ucm(:));
end
disp('The POD DEIM testing error is')
rer



%% Ratio of savings
vpods = numel(u);
tpods = numel(ulst{1}) + numel(ulst{2}) + numel(ulst{3});
disp('Percentage of storage costs is')
(tpods/vpods)*100

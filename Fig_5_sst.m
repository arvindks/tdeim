%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Arvind K.\ Saibaba 
% North Carolina State Unviersity
% If you use this code in any form, please cite "Tensor-based flow 
% reconstruction from optimally located sensor measurements" 
% by M. Farazmand and A. Saibaba, J. Fluid Mech. (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear; close all; clc
datpath = '../DATA/';
figpath = '..//figures/';
path(path, '../MATLAB')
[ Lat, Lon, time, mask, sst ] = read_data_enso( [datpath,'sst.wkmean.1990-present.nc'],...
    [datpath,'lsmask.nc'] );
rng(0);


% each element of time array is a new week, in units of days
t0 = datetime(1800,1,1,0,0,0) + days(time(1));
tfin = datetime(1800,1,1,0,0,0) + days(time(end));

[m,n,p] = size(sst);
N = m*n;
X = zeros(N,p);

M = N;
Y = zeros(M,length(time));
for i=1:length(time)
   snapshot = reshape(sst(:,:,i),N,1);
   ind = find(mask == 1);
   Y(ind,i) = snapshot(ind);

end


% train on first 16 years
Iord = 1:length(time);
Itrain = Iord(1:1200); %Iord(1:52*16);
Itest = Iord(~ismember(Iord,Itrain));

Train = Y(:,Itrain);
meansst = mean(Train,2);
Train = bsxfun(@minus,Train,meansst);




test_ind = 10;

%% Tensor based sensor placement
Ttrain = reshape(Train, 360, 180, []);
indlst = [50,25];
%figure,
for i = 1:2
   ind = indlst(i);
   Mj = tens2mat(Ttrain,i);
   [u,s,v] = svd(Mj, 'econ');
   uj = u(:,1:ind);
   
   [~,~,pj] = qr(uj',0);
   plst{i} = pj(1:ind);
   ulst{i} = uj/uj(pj(1:ind),:);
end


x = Y(:,Itest(test_ind))-meansst;
Mean = reshape(meansst, 360, 180);
XX = reshape(x,360,180);
u1 = ulst{1}; u2 = ulst{2}; p1 = plst{1}; p2 = plst{2};
XXr = mask.*(u1*XX(p1,p2)*u2');
[I,J] = find(mask == 1);
cleanind = [];
for i = 1:length(p1)
    for j = 1:length(p2)
        ind = find(I == p1(i) & J == p2(j));
        if ~isempty(ind)
            cleanind = [cleanind; p1(i), p2(j)];
        end
    end
end
xs = cleanind(:,1); ys = cleanind(:,2);




%% POD-DEIM approach
[Psi,S,V] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);

r = length(xs); % Take same number of sensors as the tensors
[~,~,pivot] = qr(Psi(:,1:r)','vector');
sensors = pivot(1:r);

x = Y(:,Itest(test_ind))-meansst;
Mean = reshape(meansst, 360, 180);

xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));


%%
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


S = [Mean(:) + XX(:), Mean(:)+XXr(:), Mean(:)+xls(:)];
bounds = [-10 + min(min(S)), max(max(S))+10];


s1 = subplot(1,3,1);
true = Mean'+ XX';
true(mask'==0) = -100;
imagesc(true)
caxis(bounds)
colormap('turbo')
title('Ground Truth', 'FontSize', 16)

s2 = subplot(1,3,2);
recon = Mean' + XXr';
recon(mask'==0)=-100;
imagesc(recon)
yticks([])
caxis(bounds)
hold on
plot(xs, ys, 'w.')
title( strcat('Tensor-DEIM, RE=', num2str(100*norm(x-XXr(:))/norm(x),'%0.2f'), '%'), 'FontSize', 16)


s3 = subplot(1,3,3);
xls = mask.*reshape(xls, 360,180);
recon = Mean' + xls';
recon(mask'==0)=-100;
imagesc(recon), colorbar
yticks([])
caxis(bounds)
ind = zeros(360*180,1);
ind(sensors) = 1;
Ind = reshape(ind, 360,180);
[xsv,ysv,~] = find(Ind);
hold on
plot(xsv, ysv, 'w.')  
title( strcat('Vector-DEIM, RE=', num2str(100*norm(x-xls(:))/norm(x),'%.2f'), '%'), 'FontSize', 16)
set(gca, 'FontSize', 16)

s1Pos = get(s1,'position');
s2Pos = get(s2,'position');
s3Pos = get(s3,'position');
s3Pos(3:4) = [s2Pos(3:4)];
set(s3,'position',s3Pos);

%print -depsc sst_plot



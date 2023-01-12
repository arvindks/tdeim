
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Arvind K.\ Saibaba 
% North Carolina State Unviersity
% If you use this code in any form, please cite "Tensor-based flow 
% reconstruction from optimally located sensor measurements" 
% by M. Farazmand and A. Saibaba, J. Fluid Mech. (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clear; close all; clc

[ Lat, Lon, time, mask, sst ] = read_data_enso('./sst.wkmean.1990-present.nc',...
    './lsmask.nc' );
rng(0);


%% Locate El Nino event in 1997
jlst = [];
for j = 1:length(time)
    date = datetime(1800,1,1,0,0,0) + days(time(j));
    if month(date) == 12 && year(date) == 1997
        jlst = [jlst; j];
    end
    
end



%% Training-test split
[m,n,p] = size(sst);
N = m*n;
X = zeros(N,p);
Y = zeros(N,length(time));
for i=1:length(time)
   snapshot = reshape(sst(:,:,i),N,1);
   ind = find(mask == 1);
   Y(ind,i) = snapshot(ind);
end

% Iord = 1:length(time);
% Itrain = Iord(1:1200); 
% Itest = Iord(~ismember(Iord,Itrain));

%% Randomly generated training data
Ntrain = round(.70*length(time)); % ~70% of data used for training
Iord = randperm(length(time));
Itrain = Iord(1:Ntrain);
Itest = Iord(~ismember(Iord,Itrain));


Train = Y(:,Itrain);
meansst = mean(Train,2);
Train = bsxfun(@minus,Train,meansst);



%% Tensor based sensor placement
Ttrain = reshape(Train, 360, 180, []);
indlst = [40,20];

for i = 1:2
   ind = indlst(i);
   Mj = tens2mat(Ttrain,i);
   [u,s,v] = svd(Mj, 'econ');
   uj = u(:,1:ind);
   
   [~,~,pj] = qr(uj',0);
   plst{i} = pj(1:ind);
   ulst{i} = uj/uj(pj(1:ind),:);
end

u1 = ulst{1}; u2 = ulst{2}; p1 = plst{1}; p2 = plst{2};
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



%% POD-DEIM 
[Psi,S,~] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);

r = length(xs); % Take same number of sensors as the tensors
[~,~,pivot] = qr(Psi(:,1:r)','vector');
sensors = pivot(1:r);
Psir = Psi(:,1:r);


%%
figure, 
fsize = 18;
rect = [0,0, 14, 6];
set(gcf, 'Units', 'inches');
set(gcf, 'OuterPosition',rect);
set(gcf, 'Position', rect);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'defaultaxesfontsize', fsize);
set(gcf, 'defaulttextfontsize', fsize);
set(0, 'DefaultAxesFontName','Times New Roman');
set(0, 'defaultTextFontName','Times New Roman');


Mean = reshape(meansst, 360, 180);
for j = 1:length(jlst)
    x = Y(:,jlst(j))-meansst;
    subplot(2,4, j)
    xls = Psir*(Psir(sensors,:)\x(sensors));
    xls = mask.*reshape(xls, 360,180);
    recon = Mean' + xls';
    imagesc(recon), colorbar
    re = norm(xls(:)-x)/norm(x);
    title(strcat('RE=', num2str(100*re, '%.2f')))
    
    subplot(2,4,j+4)
    XX = reshape(x,360,180);
    XXr = mask.*(u1*XX(p1,p2)*u2');
    recon = Mean' + XXr';
    imagesc(recon), colorbar
    re = norm(XXr(:)-x)/norm(x);
    title(strcat('RE=', num2str(100*re, '%.2f')))
    
end


%%
figure,
fsize = 18;
rect = [0,0, 12, 4];
set(gcf, 'Units', 'inches');
set(gcf, 'OuterPosition',rect);
set(gcf, 'Position', rect);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'defaultaxesfontsize', fsize);
set(gcf, 'defaulttextfontsize', fsize);
set(0, 'DefaultAxesFontName','Times New Roman');
set(0, 'defaultTextFontName','Times New Roman');
j = 1;

x = Y(:,jlst(j))-meansst;
subplot(1,3,1)
trueim = mask.*reshape(Y(:,jlst(j)), 360,180);
imagesc(trueim'), colorbar
title('True', 'FontSize', 18)

subplot(1,3,2)
xls = Psir*(Psir(sensors,:)\x(sensors));
xls = mask.*reshape(xls, 360,180);
recon = Mean' + xls';
imagesc(recon-trueim'), colorbar
re = norm(xls(:)-x)/norm(x);
%title(strcat('VectorDEIM, RE=', num2str(100*re, '%.2f'),'%'),  'FontSize', 18)
title(strcat('VectorDEIM, Diff Image'),  'FontSize', 18)

subplot(1,3, 3)
XX = reshape(x,360,180);
XXr = mask.*(u1*XX(p1,p2)*u2');
recon = Mean' + XXr';
imagesc(recon-trueim'), colorbar
re = norm(XXr(:)-x)/norm(x);
title(strcat('TensorDEIM, Diff Image'),  'FontSize', 18)
print -depsc el_nino

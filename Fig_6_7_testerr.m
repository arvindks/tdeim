clear; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Arvind K.\ Saibaba 
% North Carolina State Unviersity
% If you use this code in any form, please cite "Tensor-based flow 
% reconstruction from optimally located sensor measurements" 
% by M. Farazmand and A. Saibaba, J. Fluid Mech. (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




rng(0);

%% Read data and subtract mean
[ Lat, Lon, time, mask, sst ] = read_data_enso('./sst.wkmean.1990-present.nc',...
    'lsmask.nc' );


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

% % train on first ~70% data
% Ntrain = round(.70*length(time));
% Iord = 1:length(time);
% Itrain = Iord(1:Ntrain); 
% Itest = Iord(~ismember(Iord,Itrain));

% %% Randomly generated training data
Ntrain = round(.70*length(time)); 
Iord = randperm(length(time));
Itrain = Iord(1:Ntrain);
Itest = Iord(~ismember(Iord,Itrain));


Train = Y(:,Itrain);
meansst = mean(Train,2);
Train = bsxfun(@minus,Train,meansst);

%% Training: T-DEIM and POD-DEIM
% T-DEIM
Ttrain = reshape(Train, 360, 180, []);
Mj = tens2mat(Ttrain,1);
[Ulst{1},s1,~] = svd(Mj, 'econ');
Mj = tens2mat(Ttrain,2);
[Ulst{2},s2,~] = svd(Mj, 'econ');
figure, semilogy(diag(s1)), hold on, semilogy(diag(s2))

% POD-DEIM
[Psi,S,V] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);



%% Compare T-DEIM and V-DEIM testing errors
indlst = [10, 5; 20, 10; 30, 15; 40, 20; 50, 25; 60 30]; 
lind = length(indlst(:,1));
errt = zeros(lind, length(Itest));
errp = zeros(lind, length(Itest));
projt = zeros(lind, length(Itest));
projp = zeros(lind, length(Itest));
condt = zeros(lind, 1);
condp = zeros(lind, 1);

trainerr = zeros(lind, 2);
ns = zeros(lind,1);
plst = cell(1,2);
ulst = cell(1,2);
for j = 1:lind
    tcond = 1;
    for i = 1:2
        u = Ulst{i};
        ind = indlst(j,i);
        uj = u(:,1:ind);
   
        [~,~,pj] = qr(uj',0);
        plst{i} = pj(1:ind);
        ulst{i} = uj/uj(pj(1:ind),:);
        wlst{i} = uj;
        tcond = tcond*norm(inv(uj(pj(1:ind),:)));
    end
    u1 = ulst{1}; u2 = ulst{2}; p1 = plst{1}; p2 = plst{2};
    w1 = wlst{1}; w2 = wlst{2};
    
    
    % Find the nonmasked sensors
    [I,J] = find(mask == 1);
    cleanind = [];
    for ii = 1:length(p1)
        for jj = 1:length(p2)
            ind = find(I == p1(ii) & J == p2(jj));
            if ~isempty(ind)
                cleanind = [cleanind; p1(ii), p2(jj)];
            end
        end
    end
    
    % Tensor DEIM approach
    for k = 1:length(Itest)
        x = Y(:,Itest(k))-meansst;
        XX = reshape(x,360,180);
        XXr = mask.*(u1*XX(p1,p2)*u2');
        normx = norm(x);
        errt(j,k) = norm(x-XXr(:))/normx;
        
        XXr = mask.*(w1*(w1'*XX*w2)*w2');
        projt(j,k) = norm(x-XXr(:))/normx;
    end
    condt(j) = tcond;
    ns(j) = size(cleanind,1);

    % POD-DEIM approach
    r = size(cleanind,1); % Take same number of sensors as the tensors
    [~,~,pivot] = qr(Psi(:,1:r)','vector');
    sensors = pivot(1:r);
    Psir = Psi(:,1:r);
    for k = 1:length(Itest)
        x = Y(:,Itest(k))-meansst;
        xls = Psir*(Psir(sensors,:)\x(sensors));
        xls = mask.*reshape(xls,360,180);
        normx = norm(x);
        errp(j,k) = norm(x-xls(:))/normx;
        
        projp(j,k) = norm(x-Psir*(Psir'*x))/normx;
    end
    condp(j) = norm(inv(Psir(sensors,:)));
    
end

%% Plot

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
subplot(1,3,1)
semilogy(ns,mean(projp,2),'o-','markersize', 10,'linewidth',2)
hold on
semilogy(ns,mean(projt,2),'d-','markersize', 10,'linewidth',2)
legend('Vector-DEIM','Tensor-DEIM','fontsize',16)
set(gca,'fontsize',16)
xlabel('${n_s}$','fontsize',22,'interpreter','latex')
ylabel('Relative Error','fontsize',22,'interpreter','latex')
axis([30 1120 0.0 1.5])
title('Projection Error', 'FontSize', 18)
subplot(1,3,2)
plot(ns,condp,'o-','markersize', 10,'linewidth',2)
hold on
plot(ns,condt,'d-','markersize', 10,'linewidth',2)
set(gca,'fontsize',16)
xlabel('${n_s}$','fontsize',22,'interpreter','latex')
title('Amplification factor', 'FontSize', 18)
subplot(1,3,3)
errorbar(ns,mean(errp,2),std(errp,0,2),'o-','markersize', 10,'linewidth',2)
hold on
errorbar(ns,mean(errt,2),std(errt,0,2),'d-','markersize', 10,'linewidth',2)

set(gca,'fontsize',16)
xlabel('${n_s}$','fontsize',22,'interpreter','latex')
ylabel('Relative Error','fontsize',22,'interpreter','latex')
axis([30 1120 0.0 1.2])
title('Testing Error', 'FontSize', 18)
print -depsc testing_err_rand

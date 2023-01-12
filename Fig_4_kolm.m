%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mohammad Farazmand 
% North Carolina State Unviersity
% If you use this code in any form, please cite "Tensor-based flow 
% reconstruction from optimally located sensor measurements" 
% by M. Farazmand and A. Saibaba, J. Fluid Mech. (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; close all; clc

%% add path to tensorlab
% addpath('../tensorlab_2016-03-28')

% each element of time array is a new week, in units of days
time = 1:1000;
L1=2*pi; L2=2*pi; % domain size
n1= 128; n2=128; % n1xn2 is grid size
x1 = linspace(0,L1,n1);
x2 = linspace(0,L2,n2);
[x1,x2]=meshgrid(x1,x2);


M=n1*n2;
Y = zeros(M,length(time));
for i=1:length(time)
    rfile=['turb_w_' sprintf('%4.4i',i)];
    load(rfile);
    Y(:,i) = w(:);
end


% sequential snapshots
Ntrain = round(.75*length(time)); % 75% of data used for training
Iord = 1:length(time);
Itrain = Iord(1:Ntrain);
Itest = Iord(~ismember(Iord,Itrain));

% random snapshots
% Ntrain = round(.75*length(time)); % 75% of data used for training
% Iord = 1:length(time);
% Itrain = randi([1 length(time)], Ntrain,1);
% Itest = Iord(~ismember(Iord,Itrain));

Train = Y(:,Itrain);
meansst = mean(Train,2);
Train = bsxfun(@minus,Train,meansst);


%% Tensor based sensor placement
Ttrain = reshape(Train, n2, n1, []);

ns = [6:2:20]; % numbr of sensors
Tensor_err = nan(length(Itest),length(ns));
Vector_err = nan(length(Itest),length(ns));
for j=1:length(ns)
    j
    %% tensor-DEIM approach
    indlst = [ns(j),ns(j)]; % how many sensors
    for i = 1:2
        Mj = tens2mat(Ttrain,i);
        [u,s,~] = svd(Mj, 'econ');
        s = diag(s);
        
        ind = indlst(i);
        uj = u(:,1:ind);
        
        [~,~,pj] = qr(uj',0);
        plst{i} = pj(1:ind);
        
        ulst{i} = uj/uj(pj(1:ind),:);
    end
    
    %% vector-DEIM approach
    [Psi,S,V] = svd(Train,'econ');
    [m,n] = size(Train);
    sing = diag(S);
    
    r = ns(j)^2; % Take same number of sensors as the tensors
    [~,~,pivot] = qr(Psi(:,1:r)','vector');
    sensors = pivot(1:r);
    
    %% Computing reconstruction error
    for k=1:length(Itest)
        x = Y(:,Itest(k))-meansst;
        Mean = reshape(meansst, n2, n1);
        XX = reshape(x,n2,n1);
        u1 = ulst{1}; u2 = ulst{2}; p1 = plst{1}; p2 = plst{2};
        XXr = u1*XX(p1,p2)*u2';
        
        % tensor-DEIM error
        Tensor_err(k,j)=norm(x-XXr(:))/norm(x);
        
        % vector-DEIM error
        xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));
        Vector_err(k,j) = norm(x-xls)/norm(x);
    end
    
end

%% plotting
errorbar(ns.^2,mean(Vector_err),std(Vector_err),'o-','markersize', 10,'linewidth',2)
hold on
errorbar(ns.^2,mean(Tensor_err),std(Tensor_err),'d-','markersize', 10,'linewidth',2)
set(gca,'fontsize',16)
xlabel('$n_s$','fontsize',22,'interpreter','latex')
ylabel('Relative Error','fontsize',22,'interpreter','latex')
axis([30 410 0 1.05])
legend('Vector-DEIM','Tensor-DEIM','fontsize',16)

function [u,p] = poddeimrand(S, r)
    % 
    % Computes POD modes and the sensor locations
    % 
    % Input:
    % S       :     Snapshot data (n x k)
    % r       :     target rank (r \leq min(k, n))
    % 
    % Output:
    % u       :     Basis vectros (n x r)
    % p       :     Indices for sensor locations (1 x r)
    
    % Randomized range finder
    p = 20;
    Yr = S*randn(size(S,2), r+p);
    [Q,~] = qr(Yr,0);
    [UB,~,~] = svd(Q'*S, 'econ');
    u = Q*UB(:,1:r);
    
    % Extract sensor locations
    [~,~,p] = qr(u',0);
    p = p(1:r);
end
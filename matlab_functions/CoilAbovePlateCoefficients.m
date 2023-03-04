function [Dec, Ct, Dt] = CoilAbovePlateCoefficients(layers, kappa, omega, Cs)

assert(isscalar(kappa))

L = length(layers);
Ct = zeros(L,1);
Dt = zeros(L,1);
lambda = zeros(L,1);
bottom_z = zeros(L,1);

mu0 = 4*pi * 1e-7;
for t = 1:L
    if t == 1
        bottom_z(t) = - layers{t}.thickness;
    else
        bottom_z(t) = bottom_z(t-1) - layers{t}.thickness;
    end
    lambda(t) = sqrt(kappa^2 + 1j * omega * mu0 * layers{t}.mur * layers{t}.sigma);
end

V = zeros(2,2,L);
V(:,:,end) = eye(2);
for t = (L-1):-1:1
    mur1    = layers{t  }.mur;
    mur2    = layers{t+1}.mur;
    lambda1 = lambda(t  );
    lambda2 = lambda(t+1);
    d = -bottom_z(t);
      
    T = zeros(2,2);
    T(1,1) = exp((-lambda2 + lambda1) * d) * (1 + mur1/mur2 * lambda2/lambda1) / 2;
    T(1,2) = exp(( lambda2 + lambda1) * d) * (1 - mur1/mur2 * lambda2/lambda1) / 2;
    T(2,1) = exp((-lambda2 - lambda1) * d) * (1 - mur1/mur2 * lambda2/lambda1) / 2;
    T(2,2) = exp(( lambda2 - lambda1) * d) * (1 + mur1/mur2 * lambda2/lambda1) / 2;

    V(:,:,t) = T * V(:,:,t+1);
end

denominator = ((kappa * mur1 + lambda1) * V(1,1,1) + (kappa * mur1 - lambda1) * V(2,1,1));
Dec = ((kappa * mur1 - lambda1) * V(1,1,1) + (kappa * mur1 + lambda1) * V(2,1,1)) / denominator * Cs;
Ct(end) =                                                      (2 * kappa * mur1) / denominator * Cs;
Dt(end) = 0;
for t = 1:L-1
    Ct(t) = V(1,1,t) * Ct(end);
    Dt(t) = V(2,1,t) * Ct(end);
end
    
end

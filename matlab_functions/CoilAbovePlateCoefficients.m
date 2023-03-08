function [Dec, Ct, Dt] = CoilAbovePlateCoefficients(layers, kappa, omega, Cs)

assert(ismatrix(kappa))
assert(isscalar(omega))
assert(all(size(kappa) == size(Cs)))

L = length(layers);

Ct     = zeros([size(kappa) L]);
Dt     = zeros([size(kappa) L]);
lambda = zeros([size(kappa) L]);
bottom_z = zeros(L,1);

mu0 = 4*pi * 1e-7;
for t = 1:L
    if t == 1
        bottom_z(t) = -layers{t}.thickness;
    else
        bottom_z(t) = bottom_z(t-1) - layers{t}.thickness;
    end
    lambda(:,:,t) = sqrt(kappa.^2 + 1j * omega * mu0 * layers{t}.mur * layers{t}.sigma);
end

V11 = zeros([size(kappa) L]);
V12 = zeros([size(kappa) L]);
V21 = zeros([size(kappa) L]);
V22 = zeros([size(kappa) L]);
V11(:,:,end) = 1;  V22(:,:,end) = 1;
for t = (L-1):-1:1
    mur1    = layers{t  }.mur;
    mur2    = layers{t+1}.mur;
    lambda1 = lambda(:,:,t  );
    lambda2 = lambda(:,:,t+1);
    d = -bottom_z(t);
      
    T11 = exp((-lambda2 + lambda1) * d) .* (1 + mur1/mur2 * lambda2./lambda1) / 2;
    T12 = exp(( lambda2 + lambda1) * d) .* (1 - mur1/mur2 * lambda2./lambda1) / 2;
    T21 = exp((-lambda2 - lambda1) * d) .* (1 - mur1/mur2 * lambda2./lambda1) / 2;
    T22 = exp(( lambda2 - lambda1) * d) .* (1 + mur1/mur2 * lambda2./lambda1) / 2;

    V11(:,:,t) = T11 .* V11(:,:,t+1) + T12 .* V21(:,:,t+1);
    V12(:,:,t) = T11 .* V12(:,:,t+1) + T12 .* V22(:,:,t+1);
    V21(:,:,t) = T21 .* V11(:,:,t+1) + T22 .* V21(:,:,t+1);
    V22(:,:,t) = T21 .* V12(:,:,t+1) + T22 .* V22(:,:,t+1);
end

denominator = (kappa * mur1 + lambda1) .* V11(:,:,1) + (kappa * mur1 - lambda1) .* V21(:,:,1);
Dec = ((kappa * mur1 - lambda1) .* V11(:,:,1) + (kappa * mur1 + lambda1) .* V21(:,:,1)) ./ denominator .* Cs;
Ct(:,:,end) =                                                        (2 * kappa * mur1) ./ denominator .* Cs;
Dt(:,:,end) = 0;
for t = 1:L-1
    Ct(:,:,t) = V11(:,:,t) .* Ct(:,:,end);
    Dt(:,:,t) = V21(:,:,t) .* Ct(:,:,end);
end
    
end

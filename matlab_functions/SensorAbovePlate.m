classdef SensorAbovePlate < handle
    properties(Access=private)
        z0
        sensor_axis
        sigma
        mur
        d

        N
        Ts

        cache_z
        cache_omega
        cache_data
    end

    methods
        function obj = SensorAbovePlateCart(z0, sensor_axis, sigma, mur, d, varargin)
            assert(z0 > 0);
            assert(norm(sensor_axis) > 0);
            assert(sigma >= 0)
            assert(mur >= 1)
            assert(d > 0)

            obj.z0 = z0;
            obj.sensor_axis = sensor_axis / norm(sensor_axis);           
            obj.sigma = sigma;
            obj.mur = mur;
            obj.d = d;
            
            p = inputParser;
            addParameter(p, 'TruncateWidth', 0.2);
            addParameter(p, 'SummationTerms', 2000);
            parse(p, varargin{:});
            
            obj.N = floor(p.Results.SummationTerms / 2);
            obj.Ts = p.Results.TruncateWidth / 2 / obj.N;
            
            obj.cache_z     = [];
            obj.cache_omega = [];
            obj.cache_data  = {};
        end

        function [Ax,Ay] = A(obj, x,y,z, omega)
            z = unique(z);
            if length(z) ~= 1
                error('z must be of same value')
            end
            
            for idx = 1:length(obj.cache_z)
                if obj.cache_z(idx) == z && obj.cache_omega(idx) == omega
                    Ax = obj.cache_data{idx}.Ax(x,y);
                    Ay = obj.cache_data{idx}.Ay(x,y);
                    return
                end
            end          
            
            Ts = obj.Ts;
            N = obj.N;
            ws = 2*pi / Ts;

            z0 = obj.z0;
            sensor_axis = obj.sensor_axis;
            sigma = obj.sigma;
            mur = obj.mur;
            d = obj.d;
            
            fftAx = zeros(2*N+1);
            fftAy = zeros(2*N+1);
            for n = 1:2*N+1
            for m = 1:2*N+1
                u = (n - 1 - N) * ws / (2*N+1);
                v = (m - 1 - N) * ws / (2*N+1);
                kappa = sqrt(u^2+v^2);
                C = Cs(u,v,z0,sensor_axis,omega) * R(kappa, sigma,mur,d,omega,z);
                fftAx(n,m) = C *  1j * v;
                fftAy(n,m) = C * -1j * u;
            end
            end
            fftAx = ifftshift(fftAx);
            fftAy = ifftshift(fftAy);
            
            Ax = fft2(fftAx) / Ts^2;  Ax = fftshift(Ax);
            Ay = fft2(fftAy) / Ts^2;  Ay = fftshift(Ay);
            
            [xs,ys] = ndgrid((-N:N) * Ts, (-N:N) * Ts);
            
            obj.cache_z    (end+1) = z;
            obj.cache_omega(end+1) = omega;
            obj.cache_data {end+1} = struct;
            
            obj.cache_data{end}.Ax = griddedInterpolant(xs,ys,Ax);
            obj.cache_data{end}.Ay = griddedInterpolant(xs,ys,Ay);
            
            Ax = obj.cache_data{end}.Ax(x,y);
            Ay = obj.cache_data{end}.Ay(x,y);
        end

        function [Ex,Ey] = E(obj, x,y,z, omega)
            [Ax,Ay] = A(obj, x,y,z, omega);
            Ex = -1j * omega * Ax;
            Ey = -1j * omega * Ay;
        end
    end
end



function val = Cs(u,v, z0, t, omega)
mu0 = 4*pi * 1e-7;
kappa = sqrt(u.^2 + v.^2);
I = 1j / omega;
val = mu0 * I * exp(-kappa*z0) * (1j * u * t(1) + 1j * v * t(2) - kappa * t(3)) / 2;
end

function val = R(kappa, sigma, mur, d, omega, z)
mu0 = 4*pi * 1e-7;
lambda = sqrt(kappa.^2 + 1j*omega * mu0*mur * sigma);

c1 = lambda + kappa * mur;
c2 = lambda - kappa * mur;

num = c1 .* exp(lambda * z) + c2 .* exp(-lambda * (z+2*d));
den = c1.^2 - c2.^2 .* exp(-2*lambda*d);
val = (2 * kappa * mur) .* num ./ den;
end
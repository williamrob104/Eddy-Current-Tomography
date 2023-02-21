classdef SensorAbovePlate < handle
    properties(Access=private)
        sensor
        plate
        impl
    end

    methods
        function obj = SensorAbovePlate(z0, sensor_axis, sigma, mur, d, varargin)
            assert(z0 > 0)
            assert(numel(sensor_axis) == 3)
            assert(norm(sensor_axis) > 0)
            assert(sigma >= 0)
            assert(mur >= 1)
            assert(d > 0)

            obj.sensor = struct('z0',z0, 'axis',sensor_axis/norm(sensor_axis));
            obj.plate = struct('sigma',sigma, 'mur',mur, 'd',d);
            
            p = inputParser;
            addParameter(p, 'TruncateWidth', 0.1);
            addParameter(p, 'SummationTerms', 500);
            parse(p, varargin{:});
            
            N = floor(p.Results.SummationTerms / 2);
            delta_x = p.Results.TruncateWidth / 2 / N;
            obj.impl = struct('N',N, 'delta_x',delta_x);            
            obj.impl.cache = {};
        end

        function [Ex,Ey] = E(obj, x,y,z, omega)
            z = unique(z);
            if numel(z) ~= 1
                error('z must be of same value')
            end

            for idx = 1:length(obj.impl.cache)
                if obj.impl.cache{idx}.z == z && obj.impl.cache{idx}.omega == omega
                    Ex = obj.impl.cache{idx}.Ex(x,y);
                    Ey = obj.impl.cache{idx}.Ey(x,y);
                    return
                end
            end

            z0      = obj.sensor.z0;
            t       = obj.sensor.axis;
            sigma   = obj.plate.sigma;
            mur     = obj.plate.mur;
            d       = obj.plate.d;
            delta_x = obj.impl.delta_x;
            N       = obj.impl.N;
            delta_u = 2*pi / delta_x / (2*N+1);
            mu0     = 4*pi * 1e-7;

            [u,v] = ndgrid((-N:N) * delta_u, (-N:N) * delta_u);
            kappa = sqrt(u.^2 + v.^2);

            C = mu0/2 ./ kappa.^2 .* (u*t(1) + v*t(2) + 1j*kappa*t(3)) .* exp(-kappa*z0) .* R(kappa,sigma,mur,d,omega,z);
            fftEx = C .* -v;
            fftEy = C .*  u;
            fftEx = ifftshift(fftEx);
            fftEy = ifftshift(fftEy);

            % TODO: deal with the case omega=0, kappa=0, fftEx(1,1) neq 0
            fftEx(1,1) = 0;
            fftEy(1,1) = 0;
            
            Ex = ifft2(fftEx) / delta_x^2;  Ex = fftshift(Ex);
            Ey = ifft2(fftEy) / delta_x^2;  Ey = fftshift(Ey);
            [xs,ys] = ndgrid((-N:N) * delta_x, (-N:N) * delta_x);        
            
            cache.z     = z;
            cache.omega = omega;            
            cache.Ex    = griddedInterpolant(xs,ys,Ex);
            cache.Ey    = griddedInterpolant(xs,ys,Ey);
            obj.impl.cache{end+1} = cache;
            
            Ex = cache.Ex(x,y);
            Ey = cache.Ey(x,y);
        end
    end
end

function val = R(kappa, sigma, mur, d, omega, z)
mu0 = 4*pi * 1e-7;
lambda = sqrt(kappa.^2 + 1j*omega * mu0*mur * sigma);

c1 = lambda + kappa * mur;
c2 = lambda - kappa * mur;

mask = c1 == 0 & c2 == 0;
c1(mask) = 1 + mur;
c2(mask) = 1 - mur;

num = c1 .* exp(lambda * z) + c2 .* exp(-lambda * (z+2*d));
den = c1.^2 - c2.^2 .* exp(-2*lambda*d);
val = 2 * mur * num ./ den;
val(~mask) = val(~mask) .* kappa(~mask);
end
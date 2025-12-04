function [u2x,u2y,u2z] = secondOrderVelocity(A,phi,omega,k,theta,h,x,z,t)
    % second‐order velocity u^{(2)} (Eqs.4.19–4.24), unidirectional
    N     = numel(omega);
    u2x   = zeros(size(t));  u2y = u2x;  u2z = u2x;
    epsw  = 1e-6;
    for n = 1:N
        for m = 1:N
            % sum‐frequency
            Bp   = computeBplus(n,m,omega,k,theta,h);
            wsum = omega(n)+omega(m);
            if abs(Bp)>epsw
                xUp = Bp*(k(n)*cos(theta(n))+k(m)*cos(theta(m)));
                yUp = Bp*(k(n)*sin(theta(n))+k(m)*sin(theta(m)));
                ksum= abs(k(n)*exp(1i*theta(n))+k(m)*exp(1i*theta(m)));
                zUp = 1i*Bp*ksum*tanh(ksum*(h+z));
                phaseP = cos(wsum*t + phi(n)+phi(m));
                u2x = u2x + A(n)*A(m)* xUp .* phaseP;
                u2y = u2y + A(n)*A(m)* yUp .* phaseP;
                u2z = u2z + real(A(n)*A(m)* zUp .* phaseP);
            end
            % diff‐frequency
            Bm    = computeBminus(n,m,omega,k,theta,h);
            wdiff = omega(n)-omega(m);
            if abs(Bm)>epsw
                xUm = Bm*(k(n)*cos(theta(n))-k(m)*cos(theta(m)));
                yUm = Bm*(k(n)*sin(theta(n))-k(m)*sin(theta(m)));
                kdiff = abs(k(n)*exp(1i*theta(n)) - k(m)*exp(1i*theta(m)));
                zUm = 1i*Bm*kdiff*tanh(kdiff*(h+z));
                phaseM = cos(wdiff*t + phi(n)-phi(m));
                u2x = u2x + A(n)*A(m)* xUm .* phaseM;
                u2y = u2y + A(n)*A(m)* yUm .* phaseM;
                u2z = u2z + real(A(n)*A(m)* zUm .* phaseM);
            end
        end
    end
end
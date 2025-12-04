function [a2x,a2y,a2z] = secondOrderAcceleration(A,phi,omega,k,theta,h,x,z,t)
    % second‐order acceleration a^{(2)} via analytic d/dt of u2
    N      = numel(omega);
    a2x    = zeros(size(t));  a2y = a2x;  a2z = a2x;
    epsw   = 1e-6;
    for n = 1:N
        for m = 1:N
            % sum‐frequency
            wsum = omega(n)+omega(m);
            if abs(wsum)>epsw
                Bp   = computeBplus(n,m,omega,k,theta,h);
                xUp = Bp*(k(n)*cos(theta(n))+k(m)*cos(theta(m)));
                yUp = Bp*(k(n)*sin(theta(n))+k(m)*sin(theta(m)));
                ksum= abs(k(n)*exp(1i*theta(n))+k(m)*exp(1i*theta(m)));
                zUp = 1i*Bp*ksum*tanh(ksum*(h+z));
                phaseP = sin(wsum*t + phi(n)+phi(m));
                a2x = a2x - A(n)*A(m)* xUp .* (wsum.*phaseP);
                a2y = a2y - A(n)*A(m)* yUp .* (wsum.*phaseP);
                a2z = a2z - real(A(n)*A(m)* zUp .* (wsum.*phaseP));
            end
            % diff‐frequency
            wdiff = omega(n)-omega(m);
            if abs(wdiff)>epsw
                Bm    = computeBminus(n,m,omega,k,theta,h);
                xUm = Bm*(k(n)*cos(theta(n))-k(m)*cos(theta(m)));
                yUm = Bm*(k(n)*sin(theta(n))-k(m)*sin(theta(m)));
                kdiff= abs(k(n)*exp(1i*theta(n)) - k(m)*exp(1i*theta(m)));
                zUm = 1i*Bm*kdiff*tanh(kdiff*(h+z));
                phaseM = sin(wdiff*t + phi(n)-phi(m));
                a2x = a2x - A(n)*A(m)* xUm .* (wdiff.*phaseM);
                a2y = a2y - A(n)*A(m)* yUm .* (wdiff.*phaseM);
                a2z = a2z - real(A(n)*A(m)* zUm .* (wdiff.*phaseM));
            end
        end
    end
end
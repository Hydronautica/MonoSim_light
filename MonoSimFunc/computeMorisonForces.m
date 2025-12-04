function [Fm_wave, eta2] = computeMorisonForces(z, Nfreq, omegaVec, A_m, phiRand, k_m_vec, ...
                                                cutoff, D1, D2, Cd, Cm, rhoW, t, h, ...
                                                Hs, Tp, Le, secondOrder, irregular)
    g      = 9.81;
    nNode  = numel(z);
    Fm_wave = zeros(2*nNode, numel(t));
    % right after you unpack Nfreq, etc.
    thetaVec = zeros(1, Nfreq);   % all waves traveling in x–direction

    for i = 1:nNode
        z_b = z(i);
        if z_b < 0 || z_b > h, continue; end
        z_w = z_b - h;

        % ——— First‐order U & Udot (regular vs irregular) —————
        if ~irregular
            omega0 = 2*pi/Tp;
            k0     = omega0^2/g/tanh(omega0^2/g/tanh(omega0*h)*h);
            a      = Hs/2;
            C1     = cosh(k0*(z_w+h))/sinh(k0*h);
            U1     = a*omega0*C1.*cos(omega0*t);
            Udot1  = -a*omega0^2*C1.*sin(omega0*t);
            U2 = zeros(size(t)); Udot2 = zeros(size(t));

            if secondOrder
                % if you need 2nd-order velocity in Morison forces,
                % you could call a secondOrderVelocity() here instead.
                % For now, keep U2=0 so only η² is returned.
            end

            U    = U1 + U2;
            Udot = Udot1 + Udot2;
        else
            U    = zeros(size(t));
            Udot = zeros(size(t));
            for m = 1:Nfreq
                decay = cosh(k_m_vec(m)*(z_w + h)) / sinh(k_m_vec(m)*h);
                U    = U + A_m(m)*decay.*cos(omegaVec(m)*t + phiRand(m));
                Udot = Udot - A_m(m)*omegaVec(m)*decay.*sin(omegaVec(m)*t + phiRand(m));
            end

            % --- Add second-order irregular if requested ---
            if secondOrder
                % get the 2nd-order horizontal velocity & acceleration at x=0, z=z_b
                [u2x, ~, ~] = secondOrderVelocity( A_m, phiRand, omegaVec, k_m_vec, ...
                                                   thetaVec, h, 0, z_b, t );
                [a2x, ~, ~] = secondOrderAcceleration( A_m, phiRand, omegaVec, k_m_vec, ...
                                                       thetaVec, h, 0, z_b, t );
                U    = U    + u2x;
                Udot = Udot + a2x;
            end


        % ——— Morison drag + inertia ————————————————
        D    = (z_b <= cutoff)*D1 + (z_b > cutoff)*D2;
        A_cs = pi*(D^2)/4;
        Fd   = 0.5*rhoW*Cd*D.*abs(U).*U;
        Fi   = rhoW*Cm*A_cs.*Udot;

        Fm_wave(2*(i-1)+1, :) = (Fd + Fi).*Le;
    end

    % ——— Now get η² at z=0 using our helper —————————
    if secondOrder
        eta2 = secondOrderElevation(A_m, phiRand, omegaVec, k_m_vec, h, t);
    else
        eta2 = zeros(size(t));
    end
    end
end
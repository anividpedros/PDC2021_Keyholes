function [moid] = ComputeMoid(A,B)
    %% New orbital elements: A frame
    % incliA = 0, omegaA = 0, argpA = 0
    [incliB, omegaB, argpB] = RefFrame (A,B);

    %% Scanning orbits
    % Scanning one full revolution of meridional plane to find local minima
    % First guess for MOID is set to be big ~1e6AU

    % Angle sweeping steps
    cstep = 0.12; %rad - based on Wisnioski paper ~ 0.12rad
    stepini = 0.07; %rad - for initial tuning step
    stepfin = 1e-5; %rad - for final step of first tuning
    stepmin = 1e-14; %rad - threshold step for secondtuning

    % Initial Guess
    trueB = -2 * cstep;
    moid = 1e6;
    dist_o = 1e6;
    vdis = ones(4,1)*1e6;

    for i = 1:2
        % First triplet
        rB = B.sma * (1-B.e.^2) / (1+B.e*cos(trueB)); %compute the radius for B
        xB = rB * (cos(omegaB)*cos(argpB+trueB) - sin(omegaB)*sin(argpB+trueB)*cos(incliB));
        yB = rB * (sin(omegaB)*cos(argpB+trueB) + cos(omegaB)*sin(argpB+trueB)*cos(incliB));
        zB = rB * sin(argpB+trueB)*sin(incliB);

        rhoB = sqrt(xB.^2+yB.^2);
        L = atan2(yB,xB);

        rA = A.sma * (1-A.e.^2) / (1+A.e*cos(L)); %compute the radius for A
        rA2 = A.sma * (1-A.e.^2) / (1-A.e*cos(L));

        if abs(rhoB-rA)>abs(rhoB+rA)
            rA = rA2;
            L = L - pi;
            diff = rhoB + rA2;
        else
            diff = rhoB - rA;
        end

        D = zB.^2+(diff).^2; % square of the distance
        % storing
        if i == 1
            dist_oo = D;
        else
            dist_o = D;
            trueB_o = trueB;
            L_o = L;
        end
        trueB = trueB + cstep;
    end
    % Full revolution
    N = 0; % number of minima
    dmin = D;

    while trueB < (2*pi + cstep)
        rB = B.sma * (1-B.e.^2) / (1+B.e*cos(trueB)); %compute the radius for B
        xB = rB * (cos(omegaB)*cos(argpB+trueB) - sin(omegaB)*sin(argpB+trueB)*cos(incliB));
        yB = rB * (sin(omegaB)*cos(argpB+trueB) + cos(omegaB)*sin(argpB+trueB)*cos(incliB));
        zB = rB * sin(argpB+trueB)*sin(incliB);

        rhoB = sqrt(xB.^2+yB.^2);
        L = atan2(yB,xB);

        rA = A.sma * (1-A.e.^2) / (1+A.e*cos(L)); %compute the radius for A
        rA2 = A.sma * (1-A.e.^2) / (1-A.e*cos(L));

        if abs(rhoB-rA)>abs(rhoB+rA)
            rA = rA2;
            L = L - pi;
            diff = rhoB + rA2;
        else
            diff = rhoB - rA;
        end

        D = zB.^2+(diff).^2; % square of the distance

        if dist_o<=D && dist_o<=dist_oo
            N = N + 1;
            vtrueB(N) = trueB_o;
            vL(N) = L_o;
            vdis(N) = dist_o;
        end
        if dmin>D
            dmin = D;
        end
        dist_oo = dist_o;
        trueB_o = trueB;
        L_o = L;
        dist_o = D;
        trueB = trueB + cstep;
    end

    %% Water Procedure
    [vtrueB, vL, vdis, N] = WaterProcedure(incliB, omegaB, argpB, N, vtrueB, vL, vdis);
    %% PARALLEL TUNING
    % Move objects separately along their orbits
    % Smallest possible distance no longer meridional distance
    rBt = NaN(3,1);
    rAt = NaN(3,1);
    xBt = NaN(3,1);
    yBt = NaN(3,1);
    zBt = NaN(3,1);
    xAt = NaN(3,1);
    yAt = NaN(3,1);
    k = 1;
    while k < N+2
        if k<=N
            moid = vdis(k);
            trueB_m = vtrueB(k);
            L_m = vL(k);
            step = stepini;
            threshold = stepfin; %maybe here problem
        else
            if N == 2
                % in case of two minima are very close to each other(<1E-4 a.u.)
                % go to "water procedure"
                if (abs(vdis(1)-vdis(2))<1e-4)
                    N = 1;
                    [vtrueB, vL, vdis, N] = WaterProcedure(incliB, omegaB, argpB, N);
                    k = 1;
                else
                    if (vdis(1)<moid)
                        moid = vdis(1);
                        trueB_m = vtrueB(1);
                        L_m = vL(1);
                    end
                end
            else
                % final tuning
                for i = 1:N-1
                    if vdis(i)<moid
                        moid = vdis(i);
                        trueB_m = vtrueB(i);
                        L_m = vL(i);
                    end
                end
            end 
            step  = 2*stepini; %inital state
            threshold = stepmin; %terminal state
        end

        rBt(2) = B.sma * (1-B.e.^2) / (1+B.e*cos(trueB_m)); %compute the radius for B
        xBt(2) = rBt(2) * (cos(omegaB)*cos(argpB+trueB_m) - sin(omegaB)*sin(argpB+trueB_m)*cos(incliB));
        yBt(2) = rBt(2) * (sin(omegaB)*cos(argpB+trueB_m) + cos(omegaB)*sin(argpB+trueB_m)*cos(incliB));
        zBt(2) = rBt(2) * sin(argpB+trueB_m)*sin(incliB);

        rAt(2) = A.sma * (1-A.e.^2) / (1+A.e*cos(L_m)); %compute the radius for A
        xAt(2) = rAt(2)*cos(L_m);
        yAt(2) = rAt(2)*sin(L_m);

        aleft = true;
        aright = true;
        bleft = true;
        bright = true;
        while (step>=threshold)
            lpoints = 0;
            k1min = 1; k1max =3;
            i1min = 1; i1max =3;
            calc1 = false; calc2 = false; calc3 = false; calc4 = false;
            if (bleft)
                rBt(1) = B.sma * (1-B.e.^2) / (1+B.e*cos(trueB_m-step)); %compute the radius for B
                xBt(1) = rBt(1) * (cos(omegaB)*cos(argpB+trueB_m-step) - sin(omegaB)*sin(argpB+trueB_m-step)*cos(incliB));
                yBt(1) = rBt(1) * (sin(omegaB)*cos(argpB+trueB_m-step) + cos(omegaB)*sin(argpB+trueB_m-step)*cos(incliB));
                zBt(1) = rBt(1) * sin(argpB+trueB_m-step)*sin(incliB);
                lpoints = lpoints + 1;
            end
            if (bright)
                rBt(3) = B.sma * (1-B.e.^2) / (1+B.e*cos(trueB_m+step)); %compute the radius for B
                xBt(3) = rBt(3) * (cos(omegaB)*cos(argpB+trueB_m+step) - sin(omegaB)*sin(argpB+trueB_m+step)*cos(incliB));
                yBt(3) = rBt(3) * (sin(omegaB)*cos(argpB+trueB_m+step) + cos(omegaB)*sin(argpB+trueB_m+step)*cos(incliB));
                zBt(3) = rBt(3) * sin(argpB+trueB_m+step)*sin(incliB);
                lpoints = lpoints + 1;
            end
            if (aleft)
                rAt(1) = A.sma * (1-A.e.^2) / (1+A.e*cos(L_m-step)); %compute the radius for A
                xAt(1) = rAt(1)*cos(L_m-step);
                yAt(1) = rAt(1)*sin(L_m-step);
                lpoints = lpoints + 1;
            end
            if (aright)
                rAt(3) = A.sma * (1-A.e.^2) / (1+A.e*cos(L_m+step)); %compute the radius for A
                xAt(3) = rAt(3)*cos(L_m+step);
                yAt(3) = rAt(3)*sin(L_m+step);
                lpoints = lpoints + 1;
            end
        k1_t = 2; i1_t = 2;
        if lpoints == 1
            if (aleft) 
                i1max = 1;
            end
            if (aright) 
                i1min = 3;
            end
            if (bright) 
                k1min = 3;
            end
            if (bleft) 
                k1max = 1;
            end
        end

        if lpoints == 2
            if (aleft && bright) 
                calc1 = true;
            end
            if (aleft && bleft) 
                calc2 = true;
            end
            if (aright && bright) 
                calc3 = true;
            end
            if (aright && bleft) 
                calc4 = true;
            end
        end

        for k1 = k1min:k1max
            for i1 = i1min:i1max
                compute = true;
                if lpoints == 2
                    if i1 ~= 1
                        if (k1 ~= 3 && calc1)||(k1 ~= 1 && calc2)
                            compute = false;
                        end
                    end
                    if i1 ~= 3
                        if (k1 ~= 3 && calc3)||(k1 ~= 1 && calc4)
                            compute = false;               
                        end
                    end
                end
                if (k1 == 2 && i1 == 2)
                    compute = false;
                end
                if (compute)
                    Dx = xBt(k1)-xAt(i1);
                    Dy = yBt(k1)-yAt(i1);
                    Dz = zBt(k1);
                    dist = Dx*Dx+Dy*Dy+Dz*Dz;
                    if (dist<moid)
                        moid = dist;
                        k1_t = k1;
                        i1_t = i1;
                    end
                end
            end
        end

        if k1_t~=2||i1_t~=2
            aleft=false; aright=false; bleft=false; bright=false;
            if i1_t~=2
                if i1_t == 1
                    aleft = true;
                    L_m = L_m - step;
                    rAt(3) = rAt(2); xAt(3) = xAt(2); yAt(3) = yAt(2);
                    rAt(2) = rAt(1); xAt(2) = xAt(1); yAt(2) = yAt(1);
                else
                    aright = true;
                    L_m = L_m + step;
                    rAt(1) = rAt(2); xAt(1) = xAt(2); yAt(1) = yAt(2);
                    rAt(2) = rAt(3); xAt(2) = xAt(3); yAt(2) = yAt(3); 
                end
            end
            if k1_t~=2
                if k1_t == 1
                   bleft = true;
                    trueB_m = trueB_m - step;
                    rBt(3) = rBt(2); xBt(3) = xBt(2); yBt(3) = yBt(2); zBt(3) = zBt(2);
                    rBt(2) = rBt(1); xBt(2) = xBt(1); yBt(2) = yBt(1); zBt(2) = zBt(1);
                else
                    bright = true;
                    trueB_m = trueB_m + step;
                    rBt(1) = rBt(2); xBt(1) = xBt(2); yBt(1) = yBt(2); zBt(1) = zBt(2);
                    rBt(2) = rBt(3); xBt(2) = xBt(3); yBt(2) = yBt(3); zBt(2) = zBt(3);
                end
            end
        else
            aleft = true;
            aright = true;
            bleft = true;
            bright = true;
            step = step*0.15;
        end
        end
        if k <= N
            vdis(k) = moid;
            vtrueB(k) = trueB_m;
            vL(k) = L_m;
        end
        k = k+1;
    end
    moid = sqrt(moid);
end




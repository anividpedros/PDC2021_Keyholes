function Xdot = NBP_propagator(t,X,Input)
    
    t0 = Input.t0; % Ephemeris in seconds
    mu_vec = Input.GM_vec(:); % All gravitational parameters solar system
   
    % Asteroid State
    x = X(1:3);
    v = X(4:6);
    
    % SUN - 2BP
    ag=-mu_vec(1)/norm(x)^3*x;
    
    % PLANETS - 3BP 
    a3=zeros(3,1);
    for i=2:length(mu_vec)
        
        kep_planet = Input.IC_planet(:,i-1);
        state_planet = cspice_conics(kep_planet,t0+t);
        xthird = state_planet(1:3);
        a3=a3+mu_vec(i)*((xthird-x)/norm(xthird-x)^3-xthird/norm(xthird)^3);
        
    end

    Xdot(1:3) = v;
    Xdot(4:6) = ag+a3;

end
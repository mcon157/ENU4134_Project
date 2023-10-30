
    Psat = interp1(Tsat,Psat,Tm,'linear','extrap');
    vol_f = interp1(Tsat,vol_f,Tm,'linear','extrap');
    rho_f = 1/vol_f;
    vol_g = interp1(Tsat,vol_g,Tm,'linear','extrap');
    rho_g = 1/vol_g;
    h_f = interp1(Tsat,h_f,Tm,'linear','extrap');
    h_g = interp1(Tsat,h_g,Tm,'linear','extrap');
    mu_f = interp1(Tsat,mu_f,Tm,'linear','extrap');
    mu_g = interp1(Tsat,mu_g,Tm,'linear','extrap');
    k_f = interp1(Tsat,k_f,Tm,'linear','extrap');
    Pr_f = interp1(Tsat,Pr_f,Tm,'linear','extrap');

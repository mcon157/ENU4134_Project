Download all files to your "Matlab" folder
On step 4/15 (from slide 61-62 Module 12)
Current Issue - Tm is slightly incorrect

    %Numerical Method for Gap Conductance
    while 1
        htc_g = 5000;
        T_fo_old = 0;
        T_fo_new = T_ci(i) + q_lin/(pi*(R_ci+R_co)*htc_g);
        T_avg = ((T_fo_new+T_ci)/2)+273.15;
        k_g = 15.8*(10^-4)*(T_avg^0.79);
        gap_eff = R_ci - R_fo;
        htc_g = k_g/(gap_eff) + 5.67*(10^-8)*((T_fo_new^4) %Emissivity part cancels to 1

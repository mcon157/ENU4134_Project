clc; clear;

g = 9.81;
k_c = 15; %W/m*K
emis_f = 1;
emis_c = 1;
%Importing benchmark
open1 = fopen('benchmark_1_noSCB.txt','r');
open2 = fopen('benchmark_2_noSCB.txt','r');
open3 = fopen('benchmark_1_SCB.txt','r');
open4 = fopen('benchmark_2_SCB.txt','r');
s1 = fscanf(open1,'%f');
s2 = fscanf(open2,'%f');
s3 = fscanf(open3,'%f');
s4 = fscanf(open4,'%f');
%Importing prop tables 
proptable = readtable('proptable.txt');
Tsat = proptable{:,1};
Psat = proptable{:,2};
vol_f = proptable{:,3};
vol_g = proptable{:,4};
h_f = proptable{:,5};
h_g = proptable{:,6};
mu_f = proptable{:,7};
k_f = proptable{:,8};
Pr_f = proptable{:,9};
mu_g = proptable{:,10};

for scenario = 1:4
    s = [s1,s2,s3,s4];
    %Pulling variables from benchmark scenario
    s = s(:,scenario);
    L = s(1);
    Le = s(2);
    D = s(3);
    P =  s(4);
    Tm_inC = s(5);
    Tm_inK = Tm_inC+273.15;
    mdot =  s(6); %kg/s
    Pnom =  s(7); %MPa
    qlinmax =  s(8); %W/m
    D_ci =  s(9);
    D_fo =  s(10);
    sub = s(11);
    %Setting up z
    inlet = (-L/2)+L/400;
    outlet = L/2;
    z = linspace(inlet,outlet,400);
    
    %Pre-allocating
    h_feb = zeros(400,1);
    h_feb1 = h_f(Tm_inC);
    T_m = zeros(400,1);
    T_1co = zeros(400,1);
    T_co = zeros(400,1);
    T_ci = zeros(400,1);
    vol_l = zeros(400,1);
    vol_v = zeros(400,1);
    dP = zeros(400,1);
    T_fo = zeros(400,1);
    T_max = zeros(400,1);
    x_e = zeros(400,1);
    x = zeros(400,1);
    ZD = zeros(400,1);
    CHFR = zeros(400,1);
    %Boiling Flag
    boiling = 0;
    sboiling = 0;
    %Table Setup
    table = table(zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),'VariableNames',{'Cell','z','T_m','T_co','T_ci','T_fo','T_max','x','x_e','CHFR','dP'});
    
    for i = 1:400
        %Energy balance and T_m
        table.Cell(i) = i;
        z_i = z(i);
        table.z(i) = z_i;
        %delta z
        del_z = L/400;
        z_half = z_i - 0.5*(del_z);
        q_lin = qlinmax*cos((pi*z_i)/Le);
        q_linhalf = qlinmax*cos((pi*z_half)/Le);
        if i == 1 
            h_feb(i) = h_feb1 + (q_linhalf*(del_z))/mdot;
        else 
            h_feb(i) = h_feb(i-1) + (q_linhalf*(del_z))/mdot;
        end
        h_fsat = interp1(Psat,h_f,(Pnom*10^6),'linear','extrap');
        h_gsat = interp1(Psat,h_g,(Pnom*10^6),'linear','extrap');
        T_m(i) = interp1(h_f,Tsat,h_feb(i),'linear','extrap');
        T_sat = interp1(Psat,Tsat,(Pnom*10^6),'linear','extrap');
        if sub == 0
            %Boiling check
                if T_m(i) >= T_sat
                    boiling = 1;
                    T_m(i) = T_sat;
                    table.T_m(i) = T_m(i);
                    hfg = interp1(Tsat,h_g,T_m(i),'linear','extrap') - interp1(Tsat,h_f,T_m(i),'linear','extrap');
                    x_e(i) = (h_feb(i)-h_fsat)/(h_gsat-h_fsat);
                    x(i) = x_e(i);
                    table.x(i) = x(i);
                    table.x_e(i) = x_e(i);
                else
                    table.T_m(i) = T_m(i);
                end
        elseif sub == 1
            %SCB check
                if T_m(i) >= T_sat
                    boiling = 1;
                    T_m(i) = T_sat;
                    table.T_m(i) = T_m(i);
                else
                    table.T_m(i) = T_m(i);
                end
        end

        %T_co
        Dh = D*((4/pi)*((P/D)^2)-1);
        A = (P^2) - (pi/4)*(D^2);
        G = mdot/A;
        psi = 1.826*(P/D) - 1.0430;
        R_co = D/2;
        mu_l = interp1(Tsat,mu_f,T_m(i),'linear','extrap');
        Re = (G*Dh)/mu_l;
        Pr_l = interp1(Tsat,Pr_f,T_m(i),'linear','extrap');
        Nu = psi*(0.023*(Re^0.8)*(Pr_l^0.333));
        k_l = interp1(Tsat,k_f,T_m(i),'linear','extrap');
        htc = Nu*(k_l/Dh);
        T_co(i) = T_m(i)+(q_lin/(2*pi*R_co*htc));
        if boiling == 0
            table.T_co(i) = T_co(i);
        end
        if sub == 1
            %Sub-Cooled w/ different x and x_e
            if T_co(i) > T_sat || sboiling == 1
                boiling = 1;
                sboiling = 1;
                x_e(i) = (h_feb(i)-h_fsat)/(h_gsat-h_fsat);
                ZD(i) = z_i;
                ZD_1 = find(abs(ZD)>0,1,'first');
                x_eZD = x_e(ZD_1);
                x(i) = x_e(i)-x_eZD*exp((x_e(i)/x_eZD)-1);
                if x(i) > 0 
                    table.x(i) = x(i);
                    table.x_e(i) = x_e(i);
                end
            end
        end

        if boiling == 1
             %Schrock and Grossman
             hfg = interp1(Tsat,h_g,T_m(i),'linear','extrap') - interp1(Tsat,h_f,T_m(i),'linear','extrap');
             vol_l(i) = interp1(Tsat,vol_f,T_m(i),'linear','extrap');
             rho_l = 1/vol_l(i);
             vol_v(i) = interp1(Tsat,vol_g,T_m(i),'linear','extrap');
             rho_v = 1/vol_v(i);
             %mu_l = interp1(Tsat,mu_f,T_m(i),'linear','extrap');
             mu_v = interp1(Tsat,mu_g,T_m(i),'linear','extrap');
             Xtt = (((1-x(i))/(x(i)))^0.9)*((rho_v/rho_l)^0.5)*((mu_l/mu_v)^0.1);
             qll = q_lin/(2*pi*R_co);
             htc2 = htc*(7400*(qll/(G*hfg))+(1.11*Xtt^(-0.66)));
             T_co(i) = T_m(i)+(q_lin/(2*pi*R_co*htc2));
             table.T_co(i) = T_co(i);
        end
        %Bowring
        if boiling == 1 || x(i) > 0
            %hfg = interp1(Tsat,h_g,T_m(i),'linear','extrap') - interp1(Tsat,h_f,T_m(i),'linear','extrap');
            pr = 0.145*Pnom;
            n = 2 - 0.5*pr;
            if pr > 1 
                F1 = (pr^(-0.368))*exp(0.648*(1-pr));
                F2 = F1*(pr^(-0.448)*exp(0.245*(1-pr)))^-1;
                F3 = pr^(0.219);
                F4 = F3*(pr^1.649);
            else
                F1 = (1/1.917)*((pr^(18.942))*exp(20.89*(1-pr))+0.917);
                F2 = 1.309*F1*(pr^(1.316)*exp(2.444*(1-pr))+0.309)^-1;
                F3 = (1/1.667)*((pr^(17.023))*exp(16.658*(1-pr))+0.667);
                F4 = F3*(pr^1.649);
            end 
            A = (2.317*((hfg*Dh*G)/4)*F1)/(1+0.0143*F2*G*sqrt(Dh));
            B = (G*Dh)/4;
            C = (0.077*F3*Dh*G)/(1+(0.347*F4*(G/1356)^n));
            qcr = psi*(A - B*hfg*x(i))/C; %in W/m^2
            CHFR(i) = qcr/qll;
            if x(i) > 0
                table.CHFR(i) = CHFR(i);
            end
        end
        %Pressure Drop
        f_lo = (Re^-0.18)*(0.1339+(0.09059*((P/D)-1))+(-0.09926*((P/D)-1)^2));
        %Single-phase, no acceleration loss:
        if boiling == 0
            vol_l(i) = interp1(Tsat,vol_f,T_m(i),'linear','extrap');
            rho_l = 1/vol_l(i);
            dPg = g*rho_l;
            dP(i) = (f_lo*(1/Dh)*((G^2)/(2*rho_l))+dPg)*del_z;
            table.dP(i) = dP(i);
        %Two-phase with acceleration loss: HEM
        elseif boiling == 1
            alpha = 1/(1+((1-x(i))/x(i))*(rho_v/rho_l));
            rho_m = alpha*rho_v + (1-alpha)*rho_l;
            dPg = g*rho_m; %assuming theta = 0
            volfg = interp1(Tsat,vol_g,T_m(i),'linear','extrap') - interp1(Tsat,vol_f,T_m(i),'linear','extrap');
            del_x = abs(x(i) - x(i-1));
            dPa = (G^2)*(del_x/del_z)*volfg;
            dPf = f_lo*(1/Dh)*((G^2)/(2*rho_m));
            dP(i) = (dPf+dPa+dPg)*del_z;
            table.dP(i) = dP(i);
        end
        %T_ci
        R_ci = D_ci/2;
        T_ci(i) = T_co(i) + (q_lin/(2*pi*k_c))*log(R_co/R_ci);
        table.T_ci(i) = T_ci(i);
    
        %Numerical Method for Gap Conductance
        R_fo = D_fo/2;
        j = 1;
        htc_g = 5000;
        T_fo_o = 0;
        T_fo_n = T_ci(i) + q_lin/(pi*(R_ci+R_fo)*htc_g);
        while j == 1
            T_fo_nK = T_fo_n + 273.15;
            T_ciK = T_ci(i) + 273.15;
            T_avg = (T_fo_nK+T_ciK)/2;
            k_g = 15.8*(10^-4)*(T_avg^0.79);
            gap_eff = R_ci-R_fo;
            htc_g = k_g/gap_eff + 5.67*(10^-8)*(((T_fo_nK^4)-(T_ciK^4))/(T_fo_nK-T_ciK));%Emissivity part cancels to 1
            T_fo_n = T_ci(i) + q_lin/(pi*(R_ci+R_fo)*htc_g);
            compare = abs(T_fo_n-T_fo_o);
            if compare < 0.001
                T_fo(i) = T_fo_n;
                j = 0; 
            else
                T_fo_o = T_fo_n;
            end
        end
        table.T_fo(i) = T_fo(i);
        
        %Numerical Method for Fuel Centerline Temperature
        k_bar = 3;
        T_max_n = T_fo(i) + q_lin/(4*pi*k_bar);
        k = 1;
        while k == 1
            kdT_max = 3824*log(402.4+T_fo(i))+((6.1256*10^-11)/4)*((T_fo(i)+273)^4)+q_lin/(4*pi);
            kdT_check = 3824*log(402.4+T_max_n)+((6.1256*10^-11)/4)*((T_max_n+273)^4);
            T_max_o = T_max_n;
            T_max_n = T_max_o + (kdT_max-kdT_check)/100;
            compare = abs(kdT_check-kdT_max);
            if compare < 0.0001
                T_max(i) = T_max_n;
                k = 0;
            end
        end
        table.T_max(i) = T_max(i);
    end
    %Table Writing
    if scenario == 1
        sheet = 'PWRnoSCB';
    elseif scenario == 2
        sheet = 'BWRnoSCB';
    elseif scenario == 3
        sheet = 'PWRSCB';
    elseif scenario == 4
        sheet = 'BWRSCB';
    end
    filename1 = 'ENU_4134_Project.xlsx';
    writetable(table,filename1,'Sheet',sheet)
    clear table
    scenario = scenario+1;
end
disp('fin')

    

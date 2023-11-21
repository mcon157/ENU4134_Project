close all; clear;

%-----Constants-----
g = 9.81; %m/s^2
k_c = 15; %W/m-K
%Currently estimating fuel and cladding as black body for gap conductance
emis_f = 1; %-
emis_c = 1; %-

%-----Input-----
open1 = fopen('input.txt','r');
s = fscanf(open1,'%f');

%-----Proptable-----
proptable = readtable('proptable.txt');
Tsat = proptable{:,1}; %C
Psat = proptable{:,2}; %Pa
vol_f = proptable{:,3}; %m^3/kg
vol_g = proptable{:,4}; %m^3/kg
h_f = proptable{:,5}; %J/kg
h_g = proptable{:,6}; %J/kg
mu_f = proptable{:,7}; %kg/m-s
k_f = proptable{:,8}; %W/m-k
Pr_f = proptable{:,9}; %-
mu_g = proptable{:,10}; %kg/m-s
L = s(1); %m
Le = s(2); %m
D = s(3); %m
P =  s(4); %m
Tm_inC = s(5);
Tm_inK = Tm_inC+273.15; %K
mdot =  s(6); %kg/s
Pnom =  s(7); %MPa
qlinmax =  s(8); %W/m
D_ci =  s(9); %m
D_fo =  s(10); %m
sub = s(11); %no subcooled boiling = 0, subcooled boiling = 1
type = s(12); %PWR = 0, BWR = 1

%Setting up z -> Temperatures need to be evaluated at defined whole number cells
inlet = (-L/2)+L/400; %Inlet conditions (adding L/400 to shift z, allowing for use of Tm_inC) 
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
alert = zeros(400,1);

%-----Boiling Flag-----
boiling = 0;
sboiling = 0;

%-----Table Setup with variable names-----
table = table(zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1), ...
'VariableNames',{'Cell','z','T_m','T_co','T_ci','T_fo','T_max','x','x_e','CHFR','dP'});
    
for i = 1:400

    %Would need to consider dT/dz in LMFBR 
    %Could use Numerical Methods that are faster, more stable , and complex
    %Currently, many correlations are conservative 
    %If time dependence were considered:
    % 1. Everything that depends on (z) likely depends on (z,t)
    % 2. m(dot) becomes m(dot)(z,t)
    % 3. Velocity must be considered in Energy Balance
    % 4. Thermal capacity must be considered
    % 5. 2D iterative evaluation (r,z) becomes 3D (r,z,t)
    % 6. Limits on del z and del t
    % 7. Potential T_m_in(t)

%-----Energy balance and T_m-----
        table.Cell(i) = i;
        z_i = z(i);
        table.z(i) = z_i;
        %delta z
        del_z = L/400;
        %Specific change in z for linear heating (z=z(i)+0.5del_z)
        z_half = z_i - 0.5*(del_z); %Module 12, slide 39
        q_lin = qlinmax*cos((pi*z_i)/Le);
        %q_lin for energy balance taken between finite volumes (half-steps)
        q_linhalf = qlinmax*cos((pi*z_half)/Le);
        if i == 1 %Specific first cell with enthaply from Tm,in
            h_feb(i) = h_feb1 + (q_linhalf*(del_z))/mdot;%Module 12, eq. 61
        else 
            h_feb(i) = h_feb(i-1) + (q_linhalf*(del_z))/mdot;%Module 12, eq. 61
        end

%-----Saturation enthalpies and Temperature-----
        Pnom_Pa = Pnom*10^6;
        h_fsat = interp1(Psat,h_f,(Pnom_Pa),'linear','extrap');
        h_gsat = interp1(Psat,h_g,(Pnom_Pa),'linear','extrap');
        T_m(i) = interp1(h_f,Tsat,h_feb(i),'linear','extrap'); %Module 12, eq. 62
        T_sat = interp1(Psat,Tsat,(Pnom_Pa),'linear','extrap');

%-----Initial Boiling check-----
        if sub == 0  %SCB check
                if T_m(i) >= T_sat
                    boiling = 1;
                    T_m(i) = T_sat; %T_m becomes T_sat in boiling
                    table.T_m(i) = T_m(i);
                    hfg = interp1(Tsat,h_g,T_m(i),'linear','extrap') - interp1(Tsat,h_f,T_m(i),'linear','extrap');
                    x_e(i) = (h_feb(i)-h_fsat)/(h_gsat-h_fsat); %Module 12, eq. 74
                    %for no subcool, x_e = x
                    x(i) = x_e(i);
                    table.x(i) = x(i);
                    table.x_e(i) = x_e(i);
                else
                    table.T_m(i) = T_m(i);
                end
        elseif sub == 1
                if T_m(i) >= T_sat
                    boiling = 1;
                    T_m(i) = T_sat;
                    table.T_m(i) = T_m(i);
                else
                    table.T_m(i) = T_m(i);
                end
        end

%-----T_co-----
        %Weisman (1.3 > P\D > 1.1) 
        %Presser would probably be better to incorporate with its wider range
        %Dh = De (assuming infinite array)
        Dh = D*((4/pi)*((P/D)^2)-1); %Module 1, eq. 15
        A = (P^2) - (pi/4)*(D^2); % Flow Area
        G = mdot/A; %Mass Flux 
        psi = 1.826*(P/D) - 1.0430; %Module 1, eq. 52
        R_co = D/2;
        mu_l = interp1(Tsat,mu_f,T_m(i),'linear','extrap');
        Re = (G*Dh)/mu_l;
        Pr_l = interp1(Tsat,Pr_f,T_m(i),'linear','extrap');
        Nu = psi*(0.023*(Re^0.8)*(Pr_l^0.333)); %Module 1, eq. 51
        k_l = interp1(Tsat,k_f,T_m(i),'linear','extrap');
        htc = Nu*(k_l/Dh); %Module 1, eq. 55
        T_co(i) = T_m(i)+(q_lin/(2*pi*R_co*htc));
        if boiling == 0  %If no boiling, Tco is recorded
            table.T_co(i) = T_co(i);
        end

        if sub == 1
%-----Sub-Cooled w/ different x and x_e-----
            if T_co(i) > T_sat || sboiling == 1 
                %Added both flags for boiling and then to ensure subcooled boiling is consistent
                boiling = 1;
                sboiling = 1;
                x_e(i) = (h_feb(i)-h_fsat)/(h_gsat-h_fsat); %Module 12, eq. 74
                ZD(i) = z_i;
                ZD_1 = find(abs(ZD)>0,1,'first'); %z value for first instance of T_co > T_sat
                x_eZD = x_e(ZD_1); %x_e(ZD)
                x(i) = x_e(i)-x_eZD*exp((x_e(i)/x_eZD)-1); %Module 12, eq. 75
                if x(i) > 0 %Both qualities only write to table when x actually starts for consistency with desired outputs
                    table.x(i) = x(i);
                    table.x_e(i) = x_e(i);
                end
            end
        end
        if boiling == 1
%-----Schrock and Grossman for T_co, evaluated after T_m-----
             hfg = interp1(Tsat,h_g,T_m(i),'linear','extrap') - interp1(Tsat,h_f,T_m(i),'linear','extrap');
             vol_l(i) = interp1(Tsat,vol_f,T_m(i),'linear','extrap');
             rho_l = 1/vol_l(i);
             vol_v(i) = interp1(Tsat,vol_g,T_m(i),'linear','extrap');
             rho_v = 1/vol_v(i);
             mu_v = interp1(Tsat,mu_g,T_m(i),'linear','extrap');
             Xtt = (((1-x(i))/(x(i)))^0.9)*((rho_v/rho_l)^0.5)*((mu_l/mu_v)^0.1); %Module 9, eq. 5
             qll = q_lin/(2*pi*R_co); %Module 12, eq. 64
             htc2 = htc*(7400*(qll/(G*hfg))+(1.11*Xtt^(-0.66))); %Module 9, eq. 6 
             T_co(i) = T_m(i)+(q_lin/(2*pi*R_co*htc2)); %Module 12, eq. 65
             table.T_co(i) = T_co(i);
        end

%-----Bowring for CHF-----
        %Module 10, slide 16-17
        if boiling == 1 || x(i) > 0 %To ensure Bowring is evaluated
            pr = 0.145*Pnom; %in MPa
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

%-----CHFR alert flags, could potentially add more type CHFR alerts-----
        if type == 0
            if table.CHFR(i) < 1.3
                alert(i) = 1;
            end
        elseif type == 1
            if table.CHFR(i) < 1.9
                alert(i) = 1;
            end
        end

%-----Pressure Drop-----
        %Interior Cheng and Todreas
        f_lo = (Re^-0.18)*(0.1339+(0.09059*((P/D)-1))+(-0.09926*((P/D)-1)^2)); %Module 1, eq. 19
        %If more than a single channel estimated, would need other Cheng and Todreas subchannel Rod-Bundle correlation
%-----Single-phase, no acceleration loss: Cheng and Todreas-----
        if boiling == 0
            vol_l(i) = interp1(Tsat,vol_f,T_m(i),'linear','extrap');
            rho_l = 1/vol_l(i);
            dPg = g*rho_l; %assuming theta = 0, x = 0, Module 4, eq. 15
            dP(i) = (f_lo*(1/Dh)*((G^2)/(2*rho_l))+dPg)*del_z; %Module 1, eq. 17 (+ gravity)
            table.dP(i) = dP(i);

%-----Two-phase with acceleration loss: HEM-----
%Currently assuming f_TP = f_lo
%For added accuracy, would require SFM or other two-phase model
        elseif boiling == 1
            alpha = 1/(1+((1-x(i))/x(i))*(rho_v/rho_l)); %Module 2, eq. 60
            rho_m = alpha*rho_v + (1-alpha)*rho_l; %Module 2, eq. 52
            dPg = g*rho_m; %assuming theta = 0, Module 4, eq. 15
            volfg = interp1(Tsat,vol_g,T_m(i),'linear','extrap') - interp1(Tsat,vol_f,T_m(i),'linear','extrap');
            del_x = abs(x(i) - x(i-1));
            %estimation of dx/dz is (change in x)/(change in z)
            dPa = (G^2)*(del_x/del_z)*volfg; %Module 4, eq. 20
            dPf = f_lo*(1/Dh)*((G^2)/(2*rho_m)); %Module 4, eq. 24
            dP(i) = (dPf+dPa+dPg)*del_z; %Module 4, eq. 30 without Gas Compressibility
            table.dP(i) = dP(i);
        end

%-----T_ci-----
        R_ci = D_ci/2;
        T_ci(i) = T_co(i) + (q_lin/(2*pi*k_c))*log(R_co/R_ci); %Module 12, eq. 66
        table.T_ci(i) = T_ci(i);
    
%-----Numerical Method for Gap Conductance and T_fo-----
%Could add expressions for Gap Closure such as htc_g,closed or fitting equations for Figure 8-22, Figure 8-27
        R_fo = D_fo/2;
        j = 1;
        %-----Guess-----
        htc_g = 5000; %W * m^-2 * K^-1
        T_fo_o = 0;
        T_fo_n = T_ci(i) + q_lin/(pi*(R_ci+R_fo)*htc_g); %Module 11, eq. 41
        while j == 1
            %Require temperature in K
            T_fo_nK = T_fo_n + 273.15; %K
            T_ciK = T_ci(i) + 273.15; %K
            T_avg = (T_fo_nK+T_ciK)/2;
            k_g = 15.8*(10^-4)*(T_avg^0.79); %Module 11, 31
            gap_eff = R_ci-R_fo;
            sigmaSB = 5.67*(10^-8);
            %Currently estimating as a Black Body and Emissivity part cancels to 1
            htc_g = k_g/gap_eff + sigmaSB*(1/((1/emis_f)+(1/emis_c)-1))*(((T_fo_nK^4)-(T_ciK^4))/(T_fo_nK-T_ciK));%Module 11, eq. 36
            T_fo_n = T_ci(i) + q_lin/(pi*(R_ci+R_fo)*htc_g); %Module 11, eq. 41
            %-----Convergence Criterion----
            compare = abs(T_fo_n-T_fo_o);
            if compare < 0.001 
                T_fo(i) = T_fo_n; %if converged, writes T_fo to the table and the Numerical iteration loop breaks
                j = 0; 
            else
                T_fo_o = T_fo_n;
            end
        end
        table.T_fo(i) = T_fo(i);

%-----Numerical Method for Fuel Centerline Temperature and T_max-----
        %-----Guess-----
        k_bar = 3; %W * m^-1 * K^-1
        T_max_n = T_fo(i) + q_lin/(4*pi*k_bar); %Module 11, eq. 23
        jj = 1;
        while jj == 1
            kdT_max = 3824*log(402.4+T_fo(i))+((6.1256*10^-11)/4)*((T_fo(i)+273)^4)+q_lin/(4*pi); %Module 11, LHS of eq. 28
            kdT_check = 3824*log(402.4+T_max_n)+((6.1256*10^-11)/4)*((T_max_n+273)^4); %Module 11, RHS of eq. 28
            T_max_o = T_max_n;
            T_max_n = T_max_o + (kdT_max-kdT_check)/100;
            %-----Convergence Criterion----
            compare = abs(kdT_check-kdT_max);
            if compare < 0.0001
                T_max(i) = T_max_n;
                jj = 0;
            end
        end
        table.T_max(i) = T_max(i);
end

%-----Table Writing-----
%could add varying reactor types here, but need to specify eventually in input
if type == 0 && sub == 0
        sheet = 'PWRnoSCB';
elseif type == 1 && sub == 0
        sheet = 'BWRnoSCB';
elseif type == 0 && sub == 1
        sheet = 'PWRSCB';
elseif type == 1 && sub == 1
        sheet = 'BWRSCB';
end
filename1 = 'ENU_4134_Project.xlsx'; %Excel for submission
filename2 = 'ENU_4134_Project.txt'; %.txt for requirement / debugging
writetable(table,filename1,'Sheet',sheet)
writetable(table,filename2,'Delimiter',' ')
clear table %clearing table for next use, but directly after writing

%-----Table highlighting for CHFR flag-----
Excel = actxserver('Excel.Application'); %Loading Excel
WB = Excel.Workbooks.Open(fullfile(pwd, filename1),0,false); %Opening environment for table editing
for k = 2:401 %Need i-1 for alert, but k in cell number because shift in table writing for excel (skip first row for variable names)
        if type == 0 && sub == 0
            sheet = 'PWRnoSCB';
        elseif type == 1 && sub == 0
            sheet = 'BWRnoSCB';
        elseif type == 0 && sub == 1
            sheet = 'PWRSCB';
        elseif type == 1 && sub == 1
            sheet = 'BWRSCB';
        end
        if alert(k-1) == 1
            cell = "J"+num2str(k); %Select cells with alert
            WB.Worksheets.Item(sheet).Range(cell).Interior.ColorIndex = 3; %Color cells with alert
        end
end
%Save changes and quit Excel 
WB.Save();
WB.Close(false); 
Excel.Quit();

disp('fin')

    

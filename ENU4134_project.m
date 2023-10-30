clc; clear;

kc = 15; %W/m*K
emisf = 1;
emisc = 1;
%Importing benchmark
open1 = fopen('benchmark_1.txt','r');
benchmark1 = fscanf(open1,'%f');
open2 = fopen('benchmark_2.txt','r');
benchmark2 = fscanf(open2,'%f');
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


prompt1 = 'Which reactor type? [PWR/BWR]';
type = input(prompt1,'s');

while 1
    if strcmp("PWR",type)
        L = benchmark1(1);
        Le = benchmark1(2);
        D = benchmark1(3);
        P =  benchmark1(4);
        Tm_inC = benchmark1(5);
        Tm_inK = Tm_inC+273.15;
        mdot =  benchmark1(6); %kg/s
        Pnom =  benchmark1(7); %MPa
        qlinmax =  benchmark1(8); %W/m
        Dci =  benchmark1(9);
        Dfo =  benchmark1(10);
        prompt2 = 'Is there Subcooled Boiling? [Y/N]';
        yn = input(prompt2,'s');
        while 1
            if strcmp("Y",yn)
                fprintf("Selected Benchmark 1 with Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,Dci,Dfo)
                break
            elseif strcmp("N",yn)
                fprintf("Selected Benchmark 1 without Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,Dci,Dfo)
                break
            else
                fprintf("Please enter Y or N")
                prompt2 = '\nIs there Subcooled Boiling? [Y/N]';
                yn = input(prompt2,'s');
            end
        end
        break
    elseif strcmp("BWR",type)
        L = benchmark2(1);
        Le = benchmark2(2);
        D = benchmark2(3);
        P =  benchmark2(4);
        Tm_inC = benchmark2(5);
        Tm_inK = Tm_inC+273.15;
        mdot =  benchmark2(6); %kg/s
        Pnom =  benchmark2(7); %MPa
        qlinmax =  benchmark2(8); %W/m
        Dci =  benchmark2(9);
        Dfo =  benchmark2(10);
        prompt2 = 'Is there Subcooled Boiling? [Y/N]';
        yn = input(prompt2,'s');
        while 1
            if strcmp("Y",yn)
                fprintf("Selected Benchmark 2 with Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,Dci,Dfo)
                break
            elseif strcmp("N",yn)
                fprintf("Selected Benchmark 2 without Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,Dci,Dfo)
                break
            else
                fprintf("Please enter Y or N")
                prompt2 = '\nIs there Subcooled Boiling? [Y/N]';
                yn = input(prompt2,'s');
            end
        end
        break
    else 
        fprintf("Please enter a valid reactor type")
        prompt1 = '\nWhich reactor type? [PWR/BWR]';
        type = input(prompt1,'s');
    end
end

inlet = (-L/2)+L/400;
outlet = L/2;
z = linspace(inlet,outlet,400);
h_feb = zeros(400,1);
h_feb1 = h_f(Tm_inC);
T_m = zeros(400,1);
table = table(zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),zeros(400,1),'VariableNames',{'Cell','z','T_m','T_co','T_ci','T_fo','T_max','x','x_e','CHFR','dP (in this cell)'});

for i = 1:400
    table.Cell(i) = i;
    z_i = z(i);
    table.z(i) = z_i;
    del_z = (abs(inlet)+outlet)/400;
    %del_z = L/400;
    z_half = z_i + 0.5*(del_z);
    q_lin = qlinmax*cos((pi*z_i)/Le);
    %q_linhalf1 = qlinmax*cos((pi*z(1))/Le);
    q_linhalf = qlinmax*cos((pi*z_half)/Le);
    if i == 1 
        h_feb(i) = h_feb1 + (q_linhalf*(del_z))/mdot;
    else 
        h_feb(i) = h_feb(i-1) + (q_linhalf*(del_z))/mdot;
    end
    T_m(i) = interp1(h_f,Tsat,h_feb(i),'linear','extrap');
    table.T_m(i) = T_m(i);
end
filename = 'ENU_4134_Project.xlsx';
if strcmp("PWR",type) && strcmp("N",yn)
    sheet = 'PWRnoSCB';
elseif strcmp("PWR",type) && strcmp("Y",yn)
    sheet = 'PWRSCB';
elseif strcmp("BWR",type) && strcmp("N",yn)
    sheet = 'BWRnoSCB';
elseif strcmp("BWR",type) && strcmp("Y",yn)
    sheet = 'BWRSCB';
end
writetable(table,filename,'Sheet',sheet)
    

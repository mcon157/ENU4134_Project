clc; clear;

kc = 15; %W/m*K
emisf = 1;
emisc = 1;

format short
%Importing benchmark
open1 = fopen('benchmark_1.txt','r');
benchmark1 = fscanf(open1,'%f');
open2 = fopen('benchmark_2.txt','r');
benchmark2 = fscanf(open2,'%f');
%Importing prop tables 
proptable = readtable('proptable.txt');

Tsat = double(proptable(:,1));
Psat = double(proptable(:,2));
vol_f = double(proptable(:,3));
vol_g = proptable(:,4);
h_f = proptable(:,5);
h_g = proptable(:,6);
mu_f = proptable(:,7);
k_f = proptable(:,8);
Pr_f = proptable(:,9);
mu_g = proptable(:,10);


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
        qlinmax =  benchmark1(8); %kW/m
        Dci =  benchmark1(9);
        Dfo =  benchmark1(10);
        prompt2 = 'Is there Subcooled Boiling? [Y/N]';
        yn = input(prompt2,'s');
        while 1
            if strcmp("Y",yn)
                fprintf("Selected Benchmark 1 with Subcooled Boiling")
                break
            elseif strcmp("N",yn)
                fprintf("Selected Benchmark 1 without Subcooled Boiling")
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
        qlinmax =  benchmark2(8); %kW/m
        Dci =  benchmark2(9);
        Dfo =  benchmark2(10);
        prompt2 = 'Is there Subcooled Boiling? [Y/N]';
        yn = input(prompt2,'s');
        while 1
            if strcmp("Y",yn)
                fprintf("Selected Benchmark 2 with Subcooled Boiling")
                break
            elseif strcmp("N",yn)
                fprintf("Selected Benchmark 2 without Subcooled Boiling")
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
z = linspace(-L/2,L/2,400);
disp(z(1))
for i = 1:400
    Tm = 244;
    Psat = interp1(Tsat,Psat,Tm);
    vol_f = interp1(Tsat,vol_f,Tm);
    vol_g = proptable(:,4);
    h_f = proptable(:,5);
    h_g = proptable(:,6);
    mu_f = proptable(:,7);
    k_f = proptable(:,8);
    Pr_f = proptable(:,9);
    mu_g = proptable(:,10);
end

disp(Tsat)
disp(Psat)
disp(vol_f)

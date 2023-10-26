clc; clear;

kc = 15; %W/m*K
emisf = 1;
emisc = 1;

format short
open1 = fopen('benchmark_1.txt','r');
benchmark1 = fscanf(open1,'%f');
open2 = fopen('benchmark_2.txt','r');
benchmark2 = fscanf(open2,'%f');
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
        fprintf("Selected Benchmark 1")
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
        fprintf("Selected Benchmark 2")
        break
    else 
        fprintf("Please enter a valid reactor type")
        prompt1 = '\nWhich reactor type? [PWR/BWR]';
        type = input(prompt1,'s');
    end
end


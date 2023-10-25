clc;
clear;

%Benchmark 1 (PWR)
%Lengths all in m 
L1 = 3.66;
Le1 = 4.0;
D1 = 0.0095;
P1 =  0.0126;

Tm_in1C = 235;
Tm_in1K = Tm_in1C + 273.15;
mdot1 = 0.32; %kg/s
Pnom1 = 15.5; %MPa
qlinmax1 = 38; %kW/m
Dci1 = 0.0084;
Dfo1 = 0.0083;

%Benchmark 2 (BWR)
%Lengths all in m 
L2 = 3.05;
Le2 = 3.5;
D2 = 0.0102;
P2 =  0.013;

Tm_in2C = 265;
Tm_in2K = Tm_in2C + 273.15;
mdot2 = 0.18; %kg/s
Pnom2 = 6.89; %MPa
qlinmax2 = 32; %kW/m
Dci2 = 0.009;
Dfo2 = 0.0088;

kc = 15; %W/m*K
emisf = 1;
emisc = 1;


prompt1 = 'Which reactor type? [PWR/BWR]';
type = input(prompt1,'s');
while 1
    if typecmp("PWR",str)
        L = L1;
        Le = Le1;
        D = D1;
        P =  P1;

        Tm_inC = Tm_in1C ;
        Tm_inK = Tm_in1K;
        mdot = mdot1; %kg/s
        Pnom = Pnom1; %MPa
        qlinmax = qlinmax1; %kW/m
        Dci = Dci1;
        Dfo = Dfo1;
        fprintf("Selected Benchmark 1")
        break
    elseif typecmp("BWR",str)
        L = L2;
        Le = Le2;
        D = D2;
        P =  P2;

        Tm_inC = Tm_in2C ;
        Tm_inK = Tm_in2K;
        mdot = mdot2; %kg/s
        Pnom = Pnom2; %MPa
        qlinmax = qlinmax2; %kW/m
        Dci = Dci2;
        Dfo = Dfo2;
        fprintf("Selected Benchmark 2")
        break
    else 
        fprintf("Please enter a valid reactor type")
    end
end






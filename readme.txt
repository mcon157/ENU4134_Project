Download all files to your "Matlab" folder
On step 13/15 (from slide 61-62 Module 12)
Current Issue - CHFR and dP is off

Placeholder:


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
        D_ci =  benchmark1(9);
        D_fo =  benchmark1(10);
        prompt2 = 'Is there Subcooled Boiling? [Y/N]';
        yn = input(prompt2,'s');
        while 1
            if strcmp("Y",yn)
                fprintf("Selected Benchmark 1 with Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,D_ci,D_fo)
                sub = 1;
                break
            elseif strcmp("N",yn)
                fprintf("Selected Benchmark 1 without Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,D_ci,D_fo)
                sub = 0;
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
        D_ci =  benchmark2(9);
        D_fo =  benchmark2(10);
        prompt2 = 'Is there Subcooled Boiling? [Y/N]';
        yn = input(prompt2,'s');
        while 1
            if strcmp("Y",yn)
                fprintf("Selected Benchmark 2 with Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,D_ci,D_fo)
                sub = 1;
                break
            elseif strcmp("N",yn)
                fprintf("Selected Benchmark 2 without Subcooled Boiling")
                fprintf("\nL = %f\nLe = %f\nD = %f,\nP = %f\nTm_in(C) = %f\nmdot = %f\nPnom = %f\nqlinmax = %f\nDci = %f\nDfo = %f",L,Le,D,P,Tm_inC,mdot,Pnom,qlinmax,D_ci,D_fo)
                sub = 0;
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
filename1 = 'ENU_4134_Project.xlsx';

if strcmp("PWR",type) && strcmp("N",yn)
    sheet = 'PWRnoSCB';
elseif strcmp("PWR",type) && strcmp("Y",yn)
    sheet = 'PWRSCB';
elseif strcmp("BWR",type) && strcmp("N",yn)
    sheet = 'BWRnoSCB';
elseif strcmp("BWR",type) && strcmp("Y",yn)
    sheet = 'BWRSCB';
end
writetable(table,filename1,'Sheet',sheet)

warning off
echo off
clc

%% List of codes

pwatoolstart; % This function prints the initial conversation

demo_er=1;
while demo_er==1
    try
        sw=input(' Select demo: ');
        if sw==0
            fprintf('\n')
            return;
        elseif sw==1
            run PWAbasic
            echo off
            clc
            pwatoolstart;
        elseif sw==2
            run PWACrtPWA
            echo off
            clc
            pwatoolstart;
        elseif sw==3
            run PWACrtNon
            echo off
            clc
            pwatoolstart;
        elseif sw==4
            run PWAanalPWA;
            echo off
            clc
            pwatoolstart;
        elseif sw==5
            run PWAanalPWA_pwq;
            echo off
            clc
            pwatoolstart;
        elseif sw==6
            run PWAanalNon;
            echo off
            clc
            pwatoolstart;
        elseif sw==7
            run PWAanalNon_pwq;
            echo off
            clc
            pwatoolstart;
        elseif sw==8
            run PWAsyntPWA;
            echo off
            clc
            pwatoolstart;
        elseif sw==9
            run PWAsyntPWA_pwq;
            echo off
            clc
            pwatoolstart;
        elseif sw==10
            run PWAsyntNon;
            echo off
            clc
            pwatoolstart;
        elseif sw==11
            run PWAsyntNon_pwq;
            echo off
            clc
            pwatoolstart;
        elseif sw==12
            run ListofCommands;
            echo off
            clc
            pwatoolstart;
        else
            demo_er=1;
        end
    catch
        sw=input(' Select demo: ');
    end
    
end
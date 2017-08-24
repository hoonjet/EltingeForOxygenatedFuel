function y = computey(x,FARstoich,HCratio,OCratio)

    MWfuel = 12.011+1.008*HCratio+15.99*OCratio; %molecular weight of fuel

    y = zeros(1,4); %1: CO when rich or stoich, 2:CO2 when rich or stoich
                    %3: CO2 when lean, 4: O2 when lean
    
    if x>=FARstoich %rich or stoich mixture, CO exists, O2 does not exsit
        a = MWfuel/(138.01*x); %air
        R = 2*a+OCratio-1;
        K = 3.5; %Water-gas shift reaction
        
        A1 = K-1;
        A2 = (HCratio/2-R)*K+(R+1);
        A3 = -R;        
        b = (-A2+sqrt(A2^2-4*A1*A3))/(2*A1); %CO2
        c = R-b; 
        d = HCratio/2+b-R;         
        e = 0; 
        y(1,1) = (1-b)/(b+(1-b)+d+e+3.764*a);
        y(1,2) =  b/(b+(1-b)+d+e+3.764*a);        
        
    else %lean mixture, CO does not exist, O2 exists
        a = MWfuel/(138.01*x); %air
        b = 1; %CO2
        c = HCratio/2; % Hydrogen
        d = 0; % CO
        e = OCratio/2-HCratio/4+a-1; %O2
        y(1,3) = b/(b+(1-b)+d+e+3.764*a);
        y(1,4) = e/(b+(1-b)+d+e+3.764*a);       
        
    end

end
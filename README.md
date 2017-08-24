# EltingeForOxygenatedFuel
Generate Eltinge chart / coefficents for Oxygenated fuel
For those who want to know what the Eltinge Chart, please refer the paper
"Eltinge L. Fuel/air ratio and distribution from exhaust gas composition,SAE paper 680114"

Script is written by MATLAB 2014a

1. EltingeCoeffs.m
   Driver script that calculates Eltinge coefficents for given fuel compsition
   Input parameters: 
    - HCratio: Hydrogen/Carbon ratio of fuel
    - OCratio: Oxygen/Carbon ratio of fuel
    - % i.e) if fuel composition is CxHyOz, HCratio = y/x, OCratio = z/x
    - AFRstoich: stoichiometric air-fuel ratio
   Output result: 
    - Eltinge coefficent set EltingeCoeffs = [a b c d]
    Maldistribution factor Sx is given as
    - log(Sx) = a+b*log(CO)+c*log(O2)+d*(log(CO)*log(O2))
      Where CO,O2 is CO and O2 in % in the dry exhaust stream
      
2. getCoefficients.m 
    Compute h(x) and y(x) for given fuel-to-air ratio x
    - h(x): Exhaust air molar ratio for an element of charge
    - y(x): Concentration of compound in an element of exhaust
    
3. calculateMatrix.m
    Performing Eltinge calculation with F/A and Sx range
    - Followed the orignal scheme in the orignal paper above
    

    
    

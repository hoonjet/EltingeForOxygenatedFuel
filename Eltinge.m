
% Eltinge coefficient calcualtion script (available for oxygenated fuel as well)
% HCratio: Hydrogen/Carbon ratio of fuel
% OCratio: Oxygen/Carbon ratio of fuel
% i.e) if fuel composition is CxHyOz, HCratio = y/x, OCratio = z/x
% AFRstoich: stoichiometric air-fuel ratio

% For the details, please check the 
% "Eltinge L. Fuel/air ratio and distribution from exhaust gas composition,SAE paper 680114"

function EltingeCoeffs = Eltinge(HCratio,OCratio,AFRstoich)
    % ============= Part 1: compute h(x)and y(x) and get coefficients
    FARstoich = 1/AFRstoich;
    FAR = (0.5:0.01:1.25)'.*FARstoich; % Interpretation range
    [yCoeff,D5,D6,D7] = getCoefficients(FAR,FARstoich,HCratio,OCratio);

    % ============= Part 2: get the table as a function of FA and Sx
    FAstart = FARstoich*0.7; 
    FAend = FARstoich*1.3; 
    FAmatrix = FAstart:0.002:FAend; %Should avoid FARstoich
    SXmatrix = 0:0.0002:0.014;
    %Parameters can be adjusted 
    dataTable = calculateMatrix(FAmatrix,SXmatrix,yCoeff,FARstoich,D5,D6,D7);

    % ============ Part 3: make the fitting for given fuel 
    % fitting form log(Sx) = a+b*log(CO)+c*log(O2)+d*(log(CO)*log(O2))
    % log are common (base 10) logarithm
    nonZeroIndex = dataTable(:,3)>0 & dataTable(:,5)>0 & dataTable(:,7)>0;
    Sx = dataTable(nonZeroIndex,3);
    CO = dataTable(nonZeroIndex,5);
    O2 = dataTable(nonZeroIndex,7);
    logSx = log10(Sx);
    temp = ones(size(logSx));
    logCO = log10(CO);
    logO2 = log10(O2);
    logCombined = log10(CO).*log10(O2);
    A = [temp, logCO, logO2, logCombined];
    EltingeCoeffs = (A\logSx)'; %use least square fit

end



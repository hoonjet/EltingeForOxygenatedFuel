function [yCoeff, D5, D6, D7] = getCoefficients(FAR,FARstoich,HCratio,OCratio)

    % Compute h(x) and y(x) for given fuel-to-air ratio x
    % h(x): Exhaust air molar ratio for an element of charge
    % y(x): Concentration of compound in an element of exhaust
    % y1: CO - (rich), y2: CO2 (rich)
    % y3: CO2 (lean), y4: O2 (lean)

    h = zeros(size(FAR));
    y = zeros(size(FAR,1),4);
    stoichMinusFAR = FARstoich-FAR;
    for i = 1:length(FAR)
        h(i) = computeh(FAR(i),FARstoich,HCratio, OCratio);    
        y(i,:) = computey(FAR(i),FARstoich,HCratio, OCratio);    
    end

    h_Stoich = round(length(h)/2);
    h_lean = h(1:h_Stoich);
    x_lean = stoichMinusFAR(1:h_Stoich);
    coeffLean = polyfit(x_lean,h_lean,1);

    h_rich = h(h_Stoich:end);
    x_rich = stoichMinusFAR(h_Stoich:end);
    coeffRich = polyfit(x_rich,h_rich,1);

    D5 = (coeffLean(1,2)+coeffRich(1,2))/2;
    D6 = coeffLean(1,1);
    D7 = coeffRich(1,1)*-1;
    y_Stoich = find(FAR>FARstoich,1,'first');
    y_lean = y(1:y_Stoich-2,:);
    FAR_lean = FAR(1:y_Stoich-2);
    y3 = polyfit(FAR_lean,y_lean(:,3)*100,3);
    y4 = polyfit(FAR_lean,y_lean(:,4)*100,3);
    y_rich = y(y_Stoich:end,:);
    FAR_rich = FAR(y_Stoich:end,:);
    y1 = polyfit(FAR_rich,y_rich(:,1)*100,3);
    y2 = polyfit(FAR_rich,y_rich(:,2)*100,3);
    
    % yCoeff = coefficents C for y's 
    % y(i,x) = C(i,1)+C(i,2)*x+C(i,3)*x^2+C(i,4)*x^3
    yCoeff = [y1(end:-1:1);y2(end:-1:1);y3(end:-1:1);y4(end:-1:1)]';
end
function h = computeh(x,FARstoich,HCratio,OCratio)
% Compute h(x): number of mols of exhaust per mol of dry air in an
% element, where x is fuel-air ratio of an element, dry air

MWfuel = 12.011+1.008*HCratio+15.99*OCratio; %molecular weight of fuel

if x == FARstoich % stoichimetric
    h = (1+3.764*(1+HCratio/4-OCratio/2))/((1+HCratio/4-OCratio/2)*4.764);
elseif x<FARstoich % lean
    h = 1+(OCratio/2-HCratio/4)*138.01*x/(4.764*MWfuel);
else % rich
    a = MWfuel/(138.01*x);
    R = 2*a+OCratio-1;
    K = 3.5; %Water-gas shift reaction
  
    A1 = K-1;
    A2 = (HCratio/2-R)*K+(R+1);
    A3 = -R;
    
    b = (-A2+sqrt(A2^2-4*A1*A3))/(2*A1);
    c = R-b;    
    d = HCratio/2+b-R;
    h = (1+d+3.764*a)/(4.764*a);        
end


end
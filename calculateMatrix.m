function resultTable = calculateMatrix(FAmatrix,SXmatrix,C,STOICH,D5,D6,D7);

%Performing Eltinge calculation with F/A and Sx range
%Followed the orignal scheme in  "Eltinge L. Fuel/air ratio and distribution from exhaust gas composition,SAE paper 680114"

%FAmatrix: FA vectors (fuel-air ratio); should avoid stoich data
%SXmatrix: Sx vectors (Sx: mal-distribution vector)

% Fuel dependent factors
% C coefficient matrix from h(x) and y(x)
% STOICH : fuel/air stoichiometric ratio
% D5-D7 approximation coefficient 


CT = zeros(4,4);
CU = zeros(5,4);

[FAs, SXs] = meshgrid(FAmatrix,SXmatrix);
testVec = reshape([FAs,SXs],[],2);
resultTable = zeros(size(testVec,1),9);

    for condition = 1:size(testVec,1)

        FA = testVec(condition,1);
        SX = testVec(condition,2);
      
        for i=1:4
            CT(1,i) = C(1,i)+  C(2,i)*FA+  C(3,i)*FA^2+C(4,i)*FA^3;
            CT(2,i) = C(2,i)+2*C(3,i)*FA+3*C(4,i)*FA^2;
            CT(3,i) = C(3,i)+3*C(4,i)*FA;
            CT(4,i) = C(4,i);
        end

        D8 = D5 + D6*(STOICH-FA);
        D9 = D5 - D7*(STOICH-FA);

        if (SX<0.002) 
            U = STOICH-FA;
            B = abs(U);
            T = 99.999*U/B;
        else
            T = (STOICH-FA)/(SX+eps);
        end

        EXPT2 = (exp(-T*T/2.0))/(sqrt(2*pi));

        A = abs(T);
        S = T/A;

        SEXPF = exp(-T^2/2.0);
        TEXPF = SEXPF*(SX/(sqrt(6.28318)));

        % Performing normal probability integral 
        % exp(-x^2/2) = sigma ((-1)^n*x^(2n))/(2^n*n!)
        if A<3.5            
            BNPI = 1;
            ANPI = 1;
            N = 0;
            NF = 1.0;
            XSIGN = 1.0;
            while abs(BNPI)>0.0001
                N = N+1;
                XSIGN = -XSIGN;
                NF = NF*N;
                BNPI = XSIGN*T^(2.0*N)/(2^N*(2.0*N+1.0)*NF);
                ANPI = ANPI+BNPI;
            end
            TNPI = A*sqrt(2/pi)*ANPI;
        else
            %incase A is big (very low probability)
            TNPI = 1.0-(sqrt(2.0/pi)*SEXPF/A)*(1.0-1.0/T^2+3/T^4-15/T^6+105/T^8);
        end
        
        %CUT = Weighted total mols of dry exhaust per mol of dry air
        CUT = D5 + 0.5*(D6-D7)*(STOICH-FA)+(D6+D7)*((STOICH-FA)*(S*TNPI*0.5)+TEXPF);
        
        
        for i = 1:2
            CU(1,i) = (           CT(1,i)*D9)/CUT;
            CU(2,i) = (CT(1,i)*D7+CT(2,i)*D9)/CUT;
            CU(3,i) = (CT(2,i)*D7+CT(3,i)*D9)/CUT;
            CU(4,i) = (CT(3,i)*D7+CT(4,i)*D9)/CUT;
            CU(5,i) = (CT(4,i)*D7           )/CUT;
        end
        for i = 3:4
            CU(1,i) = (            CT(1,i)*D8)/CUT;
            CU(2,i) = (-CT(1,i)*D6+CT(2,i)*D8)/CUT;
            CU(3,i) = (-CT(2,i)*D6+CT(3,i)*D8)/CUT;
            CU(4,i) = (-CT(3,i)*D6+CT(4,i)*D8)/CUT;
            CU(5,i) = (-CT(4,i)*D6           )/CUT;
        end
        
        % Just for convinience
        SX2 = SX^2;
        SX3 = SX^3;
        SX4 = SX^4;
        T2 = T^2;
        T3 = T^3;

        CO=(CU(1,1)+CU(3,1)*SX2+3.0*CU(5,1)*SX4)*0.5*(1.0-S*TNPI)+(EXPT2)*(CU(2,1)*SX+CU(3,1)*SX2*T+CU(4,1)*SX3*(T2+2.0)+CU(5,1)*SX4*(T3+T*3.0));
        O2=(CU(1,4)+CU(3,4)*SX2+3.0*CU(5,4)*SX4)*0.5*(1.0+S*TNPI)-(EXPT2)*(CU(2,4)*SX+CU(3,4)*SX2*T+CU(4,4)*SX3*(T2+2.0)+CU(5,4)*SX4*(T3+T*3.0));
        CO2_1 = (CU(1,2)+CU(3,2)*SX2+3.0*CU(5,2)*SX4)*0.5*(1.0-S*TNPI)+(EXPT2)*(CU(2,2)*SX+CU(3,2)*SX2*T+CU(4,2)*SX3*(T2+2.0)+CU(5,2)*SX4*(T3+T*3.0));
        CO2_2 = (CU(1,3)+CU(3,3)*SX2+3.0*CU(5,3)*SX2)*0.5*(1.0+S*TNPI)-(EXPT2)*(CU(2,3)*SX+CU(3,3)*SX2*T+CU(4,3)*SX3*(T2+2.0)+CU(5,3)*SX4*(T3+T*3.0));
        CO2 = CO2_1+CO2_2;
        
        resultTable(condition,:) = [FA,STOICH,SX,T,CO,CO2,O2,EXPT2,TNPI];
    end

end



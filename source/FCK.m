function [c,k,o_ice]=FCK(T_T,P_soil,P_water,P_air,ab_a,ab_b)
%  Calculation of thermal conductivity and volumetric heat capacity
L          = 334.5;
k_soil     = 2.92;
k_air      = 0.025;
k_ice      = 2.256;
k_water    = 0.56;
c_soil     = 1.89;
c_water    = 4.2;
c_ice      = 1.94;
c_air      = 1.01;
[n,m]      = size(T_T);
T_phase    = 1;% phase change range
Tf         = -0.01;%frzeeing point
%%
for i = 1:n
    for j = 1:m
        a       = ab_a(i,j);
        b       = ab_b(i,j);
        TT      = T_T(i,j);
        o_soil  = P_soil(i,j);
        o_water = P_water(i,j);
        o_air   = P_air(i,j);
        
        %%
        if (TT >= Tf)
            o_ice(i,j) = 0;
            c(i,j)     = c_soil*o_soil+o_water*c_water+c_air*o_air;
            k(i,j)     = k_soil^o_soil* k_water^o_water*k_air^o_air;
            
        elseif (TT<-T_phase)
            o_waterl = a*(abs(Tf))^b;
            if o_waterl > o_water
                o_uwater   = 0;
                o_ice(i,j) = o_water;
            else
                o_uwater   = a*(abs(TT))^b;
                o_water    = o_waterl;
                o_ice(i,j) = o_water-o_uwater;
            end
            
            c(i,j) = c_soil*o_soil+o_uwater*c_water+c_air*o_air+o_ice(i,j)*c_ice;
            k(i,j) = k_soil^o_soil* k_water^o_uwater*k_air^o_air*k_ice^o_ice(i,j);
        else
            if TT > Tf
                TT = Tf;
            end
            o_waterl = a*(abs(Tf))^b;
            if o_waterl > o_water
                o_uwater   = 0;
                o_ice(i,j) = o_water;
            else
                o_uwater   = a*(abs(TT))^b;
                o_water    = o_waterl;
                o_ice(i,j) = o_water-o_uwater;
            end
            
            dx = -a*b*abs(TT)^(b-1);
            xb = dx*L;
            c(i,j) = o_uwater*c_water+o_ice(i,j)*c_ice+c_soil*o_soil+c_air*o_air+xb;
            k(i,j) = k_ice^o_ice(i,j)*k_soil^o_soil* k_water ^ o_uwater*k_air^o_air;
        end
    end
end
end
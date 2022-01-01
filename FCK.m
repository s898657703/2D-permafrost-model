function [c,k]=FCK(T_T,ipy,P_soil,P_water,P_air,ab_a,ab_b)
depth_peat = 0.1;%设置泥炭厚度
k_peat =0.92;
k_silt= 2.92;
k_air=0.025;
k_ice=2.256;%公式21
k_water=0.56;
c_silt=1.89;%公式27；
c_peat=1.58;
c_water=4.2;
c_ice=1.94;%公式23；冰容积热容
c_air=1.01;
T_phase=3;
[n,m] = size(T_T);
%%
for i=1:n%111行i
    for j = 1:m%97列j
        %%
        a = ab_a(i,j);
        b = ab_b(i,j);
        TT = T_T(i,j);
        o_soil = P_soil(i,j);
        o_water = P_water(i,j);
        o_air = P_air(i,j);
        
        %%      判断深度，确定土壤颗粒是否为peat或silt
        if  (ipy(i,1)<=depth_peat)
            c_soil=c_peat;%公式26；
            k_soil=k_peat;
        else
            c_soil=c_silt;%公式27；
            k_soil=k_silt;
        end
        
        %%
        if (TT>=0)
            
            c(i,j)=o_water*c_water+c_soil*o_soil+c_air*o_air; %公式17；
            k(i,j)=k_soil^o_soil* k_water^o_water*k_air^o_air;%公式16
            
        elseif (TT<-T_phase)
            
            o_ice= o_water;%THETA为土比例%公式18，求含冰量
            
            c(i,j)=o_ice*c_ice+c_soil*o_soil+c_air*o_air; %公式17；
            k(i,j)=k_ice^o_ice*k_soil^o_soil*k_air^o_air;%公式16
            
        else
            
            if TT>-0.05
                TT=-0.05;
            end
            
            o_uwater=a*(abs(TT))^b;%未冻水含量%公式15
            o_ice=a*(abs(0.001))^b - o_uwater;%THETA为土比例%公式18，求含冰量
            
            dx=-a*b*abs(TT)^(b-1);
            L=(333.2+4.995*TT+0.02987*(TT)^2);
            xb=dx*L;
            
%         if  (ipy(i,1)>=6)
%             if  xb>250%200
%                 xb=10;
%             end
%         else
%             if  xb>250%200
%                 xb=250;
%             end
%         end            
 
            c(i,j)=o_uwater*c_water+o_ice*c_ice+c_soil*o_soil+c_air*o_air+xb; %公式17；
            k(i,j)=k_ice^o_ice*k_soil^o_soil* k_water ^ o_uwater*k_air^o_air;%公式16
        end
        
        
    end
end
end
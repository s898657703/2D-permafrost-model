function [P_soil,P_water,P_air,A_t,A_f,B_t,B_f,se_1]=FPARA(ab_1,ab_2,ab_3,para_small,ipx,ipy,se)
%%      初始样品物理参数分配至每一层
ay1 =para_small(:,1);%读取参数值对应的深度
para_1(:,1)=ipy;
for i = 2:4
    para_1(:,i)=interp1(ay1,para_small(:,i),ipy,'Linear');
end
%%
ay2 =para_small(:,5);%读取参数值对应的深度
para_2(:,1)=ipy;
for i = 2:4
    para_2(:,i)=interp1(ay2,para_small(:,4+i),ipy,'Linear');
end
%%
ay3 =para_small(:,9);%读取参数值对应的深度
para_3(:,1)=ipy;
for i = 2:4
    para_3(:,i)=interp1(ay3,para_small(:,8+i),ipy,'Linear');
end


cy1 =se(:,1);%读取参数值对应的深度
se_1(:,1)=ipy;
for i = 2:3
    se_1(:,i)=interp1(cy1,se(:,i),ipy,'Linear');
end
se_1=round(se_1);



%%

by1 =ab_1(:,1);%读取参数值对应的深度
ab_1(:,1)=ipy;
for i = 2:5
    ab_1(:,i)=interp1(by1,ab_1(:,i),ipy,'Linear');
end

by2 =ab_2(:,1);%读取参数值对应的深度
ab_2(:,1)=ipy;
for i = 2:5
    ab_2(:,i)=interp1(by2,ab_2(:,i),ipy,'Linear');
end

by3 =ab_3(:,1);%读取参数值对应的深度
ab_3(:,1)=ipy;
for i = 2:5
    ab_3(:,i)=interp1(by3,ab_3(:,i),ipy,'Linear');
end

ab_f_1=ab_1(:,1:3);
ab_t_1=ab_1(:,[1 4 5]);

ab_f_2=ab_2(:,1:3);
ab_t_2=ab_2(:,[1 4 5]);

ab_f_3=ab_3(:,1:3);
ab_t_3=ab_3(:,[1 4 5]);

%%
[~,	NX] = size(ipx);
[NY,~] = size(ipy);
[~,nrr] = find(ipx==47.5);

DP_1 = (para_2-para_1)/47.5;%47.5
DP_2 = (para_3-para_2)/46;

DABF_1 = (ab_f_2-ab_f_1)/47.5;%47.5
DABF_2 = (ab_f_3-ab_f_2)/46;

DABT_1 = (ab_t_2-ab_t_1)/47.5;%47.5
DABT_2 = (ab_t_3-ab_t_2)/46;

for i = 1:NY
    
    P_soil(i,1)=para_1(i,2);
    P_water(i,1)=para_1(i,3);
    P_air(i,1)=para_1(i,4);
    
    A_t(i,1)=ab_t_1(i,2);
    A_f(i,1)=ab_f_1(i,2);
    B_t(i,1)=ab_t_1(i,3);
    B_f(i,1)=ab_f_1(i,3);
    
    for ii=2:nrr-1
        
        P_soil(i,ii)=para_1(i,2)+  DP_1(i,2) *(ipx(1,ii)-ipx(1,1));
        P_water(i,ii)=para_1(i,3)+  DP_1(i,3) *(ipx(1,ii)-ipx(1,1));
        P_air(i,ii)=para_1(i,4)+  DP_1(i,4) *(ipx(1,ii)-ipx(1,1));
        A_t(i,ii)=ab_t_1(i,2)+  DABT_1(i,2) *(ipx(1,ii)-ipx(1,1));
        A_f(i,ii)=ab_f_1(i,2)+  DABF_1(i,2) *(ipx(1,ii)-ipx(1,1));
        B_t(i,ii)=ab_t_1(i,3)+  DABT_1(i,3) *(ipx(1,ii)-ipx(1,1));
        B_f(i,ii)=ab_f_1(i,3)+  DABF_1(i,3) *(ipx(1,ii)-ipx(1,1));
        
    end
    P_soil(i,nrr)=para_2(i,2);
    P_water(i,nrr)=para_2(i,3);
    P_air(i,nrr)=para_2(i,4);
    
    A_t(i,nrr)=ab_t_2(i,2);
    A_f(i,nrr)=ab_f_2(i,2);
    B_t(i,nrr)=ab_t_2(i,3);
    B_f(i,nrr)=ab_f_2(i,3);
    
    for ii=nrr+1:NX-1
        P_soil(i,ii)=para_2(i,2)+  DP_2(i,2) *(ipx(1,ii)-ipx(1,nrr));
        P_water(i,ii)=para_2(i,3)+  DP_2(i,3) *(ipx(1,ii)-ipx(1,nrr));
        P_air(i,ii)=para_2(i,4)+  DP_2(i,4) *(ipx(1,ii)-ipx(1,nrr));
        
        A_t(i,ii)=ab_t_2(i,2)+  DABT_2(i,2) *(ipx(1,ii)-ipx(1,nrr));
        A_f(i,ii)=ab_f_2(i,2)+  DABF_2(i,2) *(ipx(1,ii)-ipx(1,nrr));
        B_t(i,ii)=ab_t_2(i,3)+  DABT_2(i,3) *(ipx(1,ii)-ipx(1,nrr));
        B_f(i,ii)=ab_f_2(i,3)+  DABF_2(i,3) *(ipx(1,ii)-ipx(1,nrr));
        
    end
    
    P_soil(i,NX)=para_3(i,2);
    P_water(i,NX)=para_3(i,3);
    P_air(i,NX)=para_3(i,4);
    
    A_t(i,NX)=ab_t_3(i,2);
    A_f(i,NX)=ab_f_3(i,2);
    B_t(i,NX)=ab_t_3(i,3);
    B_f(i,NX)=ab_f_3(i,3);
    
end

end






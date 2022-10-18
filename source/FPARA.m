function [P_soil,P_water,P_air,A_t,A_f,B_t,B_f,se_1]=FPARA(para_1,para_2,para_3,ipx,ipy,se,wide_left,wide_right,NX,NY)


[row_n,~] = size(para_1);
for i = 1 : row_n
    z_top = para_1(i,1);
    z_bot = para_1(i,2);
    
    z1 =  find(ipy>=z_top,1,'first');
    z2 =  find(ipy<=z_bot,1,'last');
    para1(z1:z2-1,:) = repmat(para_1(i,:),z2-z1,1);
    
end
para1(NY,:) = para1(NY-1,:);

[row_n,~] = size(para_2);
for i = 1 :row_n
    z_top = para_2(i,1);
    z_bot = para_2(i,2);
    
    z1 =  find (ipy>=z_top,1,'first');
    z2 =  find (ipy<=z_bot,1,'last');
    para2(z1:z2-1,:) = repmat(para_2(i,:),z2-z1,1);
    
end
para2(NY,:) = para2(NY-1,:);

[row_n,~] = size(para_3);
for i = 1 :row_n
    z_top = para_3(i,1);
    z_bot = para_3(i,2);
    
    z1 =  find (ipy>=z_top,1,'first');
    z2 =  find (ipy<=z_bot,1,'last');
    para3(z1:z2-1,:) = repmat(para_3(i,:),z2-z1,1);
    
end
para3(NY,:) = para3(NY-1,:);
%%      定义冻融日期
cy1 =se(:,1);%读取参数值对应的深度
se_1(:,1)=ipy;
for i = 2:3
    se_1(:,i)=interp1(cy1,se(:,i),ipy,'Linear');
end
se_1=round(se_1);


%%  含水量
ab_f_1=para1(:,6:7);
ab_t_1=para1(:,8:9);

ab_f_2=para2(:,6:7);
ab_t_2=para2(:,8:9);

ab_f_3=para3(:,6:7);
ab_t_3=para3(:,8:9);

n1 =  find(ipx>=wide_left,1,'first');
n2 =  find(ipx>=wide_right,1,'first');


for i = 1:NY
    for ii=1:n1
        A_t(i,ii)=ab_t_1(i,1);
        B_t(i,ii)=ab_t_1(i,2);
        A_f(i,ii)=ab_f_1(i,1);
        B_f(i,ii)=ab_f_1(i,2);
        P_soil(i,ii) = para1(i,3);
        P_water(i,ii) = para1(i,4);
        P_air(i,ii) = para1(i,5);
         
    end
    for ii=n1+1:n2
        
        A_t(i,ii)=ab_t_2(i,1);
        A_f(i,ii)=ab_f_2(i,1);
        B_t(i,ii)=ab_t_2(i,2);
        B_f(i,ii)=ab_f_2(i,2);
        P_soil(i,ii) = para2(i,3);
        P_water(i,ii) = para2(i,4);
        P_air(i,ii) = para2(i,5);
    end
    for ii=n2+1:NX
        
        A_t(i,ii)=ab_t_3(i,1);
        A_f(i,ii)=ab_f_3(i,1);
        B_t(i,ii)=ab_t_3(i,2);
        B_f(i,ii)=ab_f_3(i,2);
        P_soil(i,ii) = para3(i,3);
        P_water(i,ii) = para3(i,4);
        P_air(i,ii) = para3(i,5);
    end    
end
end






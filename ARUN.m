%%      验证
clc
clear
Tr1=0.0651 ; %warming rate 

filepath1 = 'E:\mat\2w_true\2D\';%Result storage path

load('T_IC.mat')%initial condition 
load('Apara.mat')%Result storage path
load('GST.mat')%uper boundary conditions 

[ipx,ipy]=Fipxy;%Gridding 

[NX,NY,T_T]=FINI(T_1,T_3,T_2,ipx,ipy);%initialization 

%%     
[NUM_node,NUM_tri,NUM_upnode,NUM_bontri,NODE] = FCALCULATOR(NX,NY);%Geometric calculation 

%% 
[XYN] = FXYN(NY,NX,NUM_node,ipx,ipy);

%%      
[IJM] = FIJM(NY,NX);

%%      
[NN,NUM_max,NNN]=FNN(IJM,NUM_node,NUM_tri);
%%      不用考虑冻土宽度，参数根据钻孔的值进行插值
[P_soil,P_water,P_air,A_t,A_f,B_t,B_f,se_1]=FPARA(ab_1,ab_2,ab_3,para_small,ipx,ipy,se);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      开始计算
tic;
h=waitbar(0,'优化计算中，请稍候！');
year_start=2018;
year_end=2030;

T_rate=0;%为升温速率
%%      利用初始温度值建立温度场
[TTT] = FTTT(T_T,NX,NY);
day_start=1;
day_end=365;
for year_run = year_start:year_end
       
  
      GST = GST + (year_run-2021)*T_rate ; 
    for  day_run = day_start:day_end

        %%      判断冻结期和融化期
        
        for ssee =1:NY
            
            if (day_run>se_1(ssee,2)&&day_run<se_1(ssee,3))
                
                ab_a(ssee,:)=A_t(ssee,:);
                ab_b(ssee,:)=B_t(ssee,:);
                
            else%freezing
                ab_a(ssee,:)=A_f(ssee,:);
                ab_b(ssee,:)=B_f(ssee,:);
            end
            
            
        end
        
        [c,k]=FCK(T_T,ipy,P_soil,P_water,P_air,ab_a,ab_b);
        
        %%      边界热流,分配不同上表面温度
        [NODE_up,Q2]=FBOUNDRY(T_T,ipx,ipy,GST,day_run,k);

        
        %%      预分配内存
        TKK = zeros(NUM_max,1);
        TNN = zeros(NUM_max,1);
        TPP = zeros(NUM_node,1);
        TSS = zeros(NUM_node,1);
        %%      开始分析每个三角形元素
        for NUM_tri_run = 1:NUM_tri%三角形单元数
            %%      求三角形单元的导热系数和容积热容
            I1 = IJM(NUM_tri_run,1);%I1、J1、M1每个三角形的三个节点号
            J1 = IJM(NUM_tri_run,2);
            M1 = IJM(NUM_tri_run,3);
            %判定节点号I1，J1，M1所在的行列号以确定其导热系数和容积热容
            [n_rowI,n_colI] = find(NODE==I1);
            [n_rowJ,n_colJ] = find(NODE==J1);
            [n_rowM,n_colM] = find(NODE==M1);
            CCI = c(n_rowI,n_colI);
            RKKI= k(n_rowI,n_colI);
            
            CCJ = c(n_rowJ,n_colI);
            RKKJ = k(n_rowJ,n_colI);
            
            CCM = c(n_rowM,n_colJ);
            RKKM = k(n_rowM,n_colJ);
            PCC = (CCI+CCJ+CCM)/3*1000000;%  求三个节点的容积热容量和导热系数的平均值
            RK = (RKKI+RKKJ+RKKM)/3*3600.0 * 24;%作为三角形单元的容积热容量和导热系数
            %%      变分计算    平面法（不对称）
            %求三角形单元的面积
            BI = XYN(J1,2)-XYN(M1,2);%插值函数的参数
            BJ = XYN(M1,2)-XYN(I1,2);%插值函数的参数
            BM = XYN(I1,2)-XYN(J1,2);%插值函数的参数
            CI = XYN(M1,1)-XYN(J1,1);%插值函数的参数
            CJ = XYN(I1,1)-XYN(M1,1);%插值函数的参数
            CM = XYN(J1,1)-XYN(I1,1);%插值函数的参数
            SS = (BI*CJ-BJ*CI)/2;%三角形单元面积
            SI = sqrt(BI^2+CI^2);%为两节点的距离
            
            FA = RK/4/SS;%对应于fi
            %求节点的导热系数
            RK11 = FA*(BI^2+CI^2);
            RK22 = FA*(BJ^2+CJ^2);
            RK33 = FA*(BM^2+CM^2);
            RK12 = FA*(BI*BJ+CI*CJ);
            RK21 = RK12;
            RK13 = FA*(BI*BM+CI*CM);
            RK31 = RK13;
            RK23 = FA*(BJ*BM+CJ*CM);
            RK32 = RK23;
            %求节点的形状矩阵
            RN11=SS*PCC/6;
            RN22=SS*PCC/6;
            RN33=SS*PCC/6;
            RN12=SS*PCC/12;
            RN21=SS*PCC/12;
            RN13=SS*PCC/12;
            RN31=SS*PCC/12;
            RN23=SS*PCC/12;
            RN32=SS*PCC/12;
            %%      判断是否为边界条件，如果是边界条件执行不同算法；
            %源程序有判断是否为边界，执行不同算法，主要是Q2,而如果直接给q2定义值，不是边界就是0 ；
            %就不用判断是否为边界条件
            PP1 = 0;
            PP2 = -Q2(NUM_tri_run,1)*SI/2;%84页，由于qv内部热流量为0，所以省略前面的内容
            PP3 = -Q2(NUM_tri_run,1)*SI/2;
            
            %%
            %前面定义TKK为0矩阵，所以在此要根据前面计算的K值重新复制到TKK矩阵
            %根据每个系数的地址确定主对角元素Kll的导热系数；
            %NN(I1)即文章中的主对角元素地址N[l]136页
            TKK(NN(I1)) = TKK(NN(I1))+RK11;
            TKK(NN(J1)) = TKK(NN(J1))+RK22;
            TKK(NN(M1)) = TKK(NN(M1))+RK33;
            TNN(NN(I1)) = TNN(NN(I1))+RN11;
            TNN(NN(J1)) = TNN(NN(J1))+RN22;
            TNN(NN(M1)) = TNN(NN(M1))+RN33;
            %%
            %NN(I1)-I1+J1为计算任意非零元素Klp的地址138页；
            %确定非0元素Klp的导热系数139页
            if (I1>J1)
                TKK(NN(I1)-I1+J1) = TKK(NN(I1)-I1+J1)+RK12;
                TNN(NN(I1)-I1+J1) = TNN(NN(I1)-I1+J1)+RN12;
            else
                TKK(NN(J1)-J1+I1) = TKK(NN(J1)-J1+I1)+RK21;
                TNN(NN(J1)-J1+I1) = TNN(NN(J1)-J1+I1)+RN21;
            end
            if (I1>M1)
                TKK(NN(I1)-I1+M1) = TKK(NN(I1)-I1+M1)+RK13;
                TNN(NN(I1)-I1+M1) = TNN(NN(I1)-I1+M1)+RN13;
            else
                TKK(NN(M1)-M1+I1) = TKK(NN(M1)-M1+I1)+RK31;
                TNN(NN(M1)-M1+I1) = TNN(NN(M1)-M1+I1)+RN31;
            end
            if (J1>M1)
                TKK(NN(J1)-J1+M1) = TKK(NN(J1)-J1+M1)+RK23;
                TNN(NN(J1)-J1+M1) = TNN(NN(J1)-J1+M1)+RN23;
            else
                TKK(NN(M1)-M1+J1) = TKK(NN(M1)-M1+J1)+RK32;
                TNN(NN(M1)-M1+J1) = TNN(NN(M1)-M1+J1)+RN32;
            end
            %%  确定热容
            TPP(I1) = TPP(I1)+PP1;
            TPP(J1) = TPP(J1)+PP2;
            TPP(M1) = TPP(M1)+PP3;
            
            
        end
        %%
        [TKK,TNN,TPP]=FCRANK(NUM_node,TPP,NUM_max,TKK,TNN);
        % %CRANK-NICOLSON
        
        %%
        %MULTI针对系数矩阵N,判断矩阵系数是否为0，是否要执行消去运算140页
        [TPP,TSS]=FMULTI(NUM_node,NNN,NN,TSS,TPP,TTT,TNN);
        %%
        for i = 1:NUM_upnode  %m19--IB%
            MRBPT0 = NODE_up(i,1);%存储的是上边界节点号
            TKK(NN(MRBPT0)) = 1.0E35;
            TPP(MRBPT0) = NODE_up(i,2)*1.0E35;
        end
        %%
        %GAUSSIAN  METHOD高斯迭代法124页
        [TTT]=FGAUSSIAN(NUM_node,NNN,NN,TSS,TPP,TKK);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        T_NXNY=zeros(NY,NX);%计算的温度存储为NY行NX列
        for i = 1:NY
            for j = 1:NX
                T_NXNY(i,j) =  TTT(NX*(NY-i)+j,1);%计算的温度存储为NY行NX列
            end
        end
        
        T_T = T_NXNY;%
        clear T_NXNY
        AJG{day_run,1}=T_T;
        str=['运行中...',num2str(year_run),'年'];
        waitbar(day_run/365,h,str)
        clear c k
    end
    year1 = sprintf('%04d',year_run);
    na = [filepath1,num2str(year1),'a.mat'];
    save(na,'AJG');
    clear AJG
end
%

close(h);
toc;
%%
%  ERR
%  drawerr


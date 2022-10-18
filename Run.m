% -------------------------------------------------------------------------
% 2D permarfost model
% main script for running the model
%
% Developed by: Wen sun 2022
%
% -------------------------------------------------------------------------

%%   
clc
clear
datapath = 'data\';
outpath  = 'out\';
%%  import soil parameter
load([datapath,'PARA.mat'])
%para_1, para_2, para_3 is three boreholes soil parameter 
% column 1:depth of layer top(unit:m)
% column 2:depth of layer bottom(unit:m)
% column 3:soil particle content (unit:m3/m3)
% column 4:water content (unit:m3/m3)
% column 5:air content (unit:m3/m3)
% column 6:fitting constant 'a' of unfrozen water in freezing state
% column 7:fitting constant 'b' of unfrozen water in freezing state
% column 8:fitting constant 'a' of unfrozen water in thawing state
% column 9:fitting constant 'b' of unfrozen water in thawing state
% se is freezing and thawing time
% column 1:depth of layer (unit:m)
% column 2:time of thawing(unit:m)
% column 3:time of freezing(unit:m)

%%  import observed soil temperature
load([datapath,'T_ini.mat'])
% T_1, T_2, T_3 is the soil temperature observed in three boreholes
% column 1:depth of layer(unit:m)
% column 2:soil temperature(unit:℃)
%%  import forcing data (GST)
load([datapath,'GST.mat'])
% GST_data is ground surface temperature
% column 1:date
% column 2-4:ground surface temperature at S1, S2, S3 site(unit:℃)
RDATE  = datestr(GST_data.date); % date
GST  = [GST_data.S1 GST_data.S2 GST_data.S3]; % ground surface temperature for three boreholes
%%  define input parameters
%                 left  right  dx
Mesh.ipx = [ 0    93.5   0.50];   
% column 1:distance of column left(unit:m)
% column 2:distance of column right(unit:m)
% column 3:horizontal mesh size between the layers(unit:m)

%                 top   bot  dy
Mesh.ipy = [ 0.1   2    0.10;...   
                  2     20   0.25;...   
                  20    30   0.5 ;...    
                  30    40   1   ; ...
                  40    60   2   ;...
                  60    100  5   ];   
% column 1:depth of layer top(unit:m)
% column 2:depth of layer bottom(unit:m)
% column 3:vertical mesh size between the layers(unit:m)
%%
begDate       = datetime(2015,10,21);  % start date for run
endDate       = datetime(2021,12,13); % end date for run
year_start    = 2015;
year_end      = 2021;
heat_flux     = 0.06;%Geothermal heat flux at the bottom
wide_left     = 22;%middle position between left two boreholes(unit:m) 
wide_right    = 70;%middle position between right two boreholes(unit:m) 
time_interval = 1;% interval [days] in which output files are written
TIME_s        = 3600.0 * 24 * time_interval;
NTB           = strmatch(datestr(begDate), RDATE); % forcing index of begDate
%%  Discretization
% Generate mesh
[ipx,ipy] = FMESH(Mesh);
NX = size(ipx,1);
NY = size(ipy,1);
%%  Analysis domain discretization
[NUM_node,NUM_tri,NUM_upnode,NUM_bontri,NODE] = FCALCULATOR(NX,NY);

%%  Assign a set of node coordinates to each node
[XYN] = FXYN(NY,NX,NUM_node,ipx,ipy);

%%  Assign a number to each triangle. Assign node numbers to the three vertices of each triangle element.
[IJM] = FIJM(NY,NX);

%%  Determine the address and maximum storage quantity of each main diagonal element of coefficient matrix in vector memory
[NN,NUM_max,NNN] = FNN(IJM,NUM_node,NUM_tri);

%%  Assign soil parameters to each layer of the analysis domain
[P_soil,P_water,P_air,A_t,A_f,B_t,B_f,se_1]=FPARA(para_1,para_2,para_3,ipx,ipy,se,wide_left,wide_right,NX,NY);

%%  initialization temperature filed
[T_T] = FINI(T_1,T_2,T_3,ipx,ipy);
[TTT] = FTTT(T_T,NX,NY);

%% run modle
tic;
h=waitbar(0,'Calculation in progress, please wait！');

day_total = NTB;
for year_run = year_start:year_end
    
    if(( rem(year_run,100)~= 0 && rem(year_run,4) == 0 )|| (rem(year_run,100) == 0 && rem(year_run,400) == 0))
        day_end = 366;
    else
        day_end = 365;
    end
    
    if year_run == 2015
        day_start = NTB;
    else
        day_start = 1;
    end
    
    for  day_run = day_start:day_end
        
%%  Determination of freezing and thawing periods
        for ssee = 1:NY
            if (day_run>se_1(ssee,2)&&day_run<se_1(ssee,3))
                ab_a(ssee,:) = A_t(ssee,:);
                ab_b(ssee,:) = B_t(ssee,:);
            else%freezing
                ab_a(ssee,:) = A_f(ssee,:);
                ab_b(ssee,:) = B_f(ssee,:);
            end
        end
        
%%  Calculation of thermal conductivity and volumetric heat capacity
        [c,k,o_ice] = FCK(T_T,P_soil,P_water,P_air,ab_a,ab_b);

%%  Determine the number and type of boundary triangle elements and set corresponding boundary conditions
        GST_day = GST(day_total,:);
        [NODE_up,Q2]=FBOUNDRY(ipx,GST_day,TIME_s,wide_left,wide_right,NX,NY,heat_flux);
%%  Calculation by triangle element
        TKK = zeros(NUM_max,1);
        TNN = zeros(NUM_max,1);
        TPP = zeros(NUM_node,1);
        TSS = zeros(NUM_node,1);
        for NUM_tri_run = 1:NUM_tri

            I1 = IJM(NUM_tri_run,1);
            J1 = IJM(NUM_tri_run,2);
            M1 = IJM(NUM_tri_run,3);
            
            [n_rowI,n_colI] = find(NODE==I1);
            [n_rowJ,n_colJ] = find(NODE==J1);
            [n_rowM,n_colM] = find(NODE==M1);
            CCI = c(n_rowI,n_colI);
            RKKI= k(n_rowI,n_colI);
            
            CCJ = c(n_rowJ,n_colI);
            RKKJ = k(n_rowJ,n_colI);
            
            CCM = c(n_rowM,n_colJ);
            RKKM = k(n_rowM,n_colJ);
            PCC = (CCI+CCJ+CCM)/3*1000000;
            RK = (RKKI+RKKJ+RKKM)/3*TIME_s;
%%  variational method     
            BI = XYN(J1,2)-XYN(M1,2);
            BJ = XYN(M1,2)-XYN(I1,2);
            BM = XYN(I1,2)-XYN(J1,2);
            CI = XYN(M1,1)-XYN(J1,1);
            CJ = XYN(I1,1)-XYN(M1,1);
            CM = XYN(J1,1)-XYN(I1,1);
            SS = (BI*CJ-BJ*CI)/2;
            SI = sqrt(BI^2+CI^2);
            FA = RK/4/SS;

            RK11 = FA*(BI^2+CI^2);
            RK22 = FA*(BJ^2+CJ^2);
            RK33 = FA*(BM^2+CM^2);
            RK12 = FA*(BI*BJ+CI*CJ);
            RK21 = RK12;
            RK13 = FA*(BI*BM+CI*CM);
            RK31 = RK13;
            RK23 = FA*(BJ*BM+CJ*CM);
            RK32 = RK23;

            RN11=SS*PCC/6;
            RN22=SS*PCC/6;
            RN33=SS*PCC/6;
            RN12=SS*PCC/12;
            RN21=SS*PCC/12;
            RN13=SS*PCC/12;
            RN31=SS*PCC/12;
            RN23=SS*PCC/12;
            RN32=SS*PCC/12;
   
            PP1 = 0;
            PP2 = -Q2(NUM_tri_run,1)*SI/2;
            PP3 = -Q2(NUM_tri_run,1)*SI/2;
            
%%  steady-state temperature filed coefficeient matrix
            TKK(NN(I1)) = TKK(NN(I1))+RK11;
            TKK(NN(J1)) = TKK(NN(J1))+RK22;
            TKK(NN(M1)) = TKK(NN(M1))+RK33;
            TNN(NN(I1)) = TNN(NN(I1))+RN11;
            TNN(NN(J1)) = TNN(NN(J1))+RN22;
            TNN(NN(M1)) = TNN(NN(M1))+RN33;

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
  
            TPP(I1) = TPP(I1)+PP1;
            TPP(J1) = TPP(J1)+PP2;
            TPP(M1) = TPP(M1)+PP3;
        end
%%  CRANK-NICOLSON
        [TKK,TNN,TPP]=FCRANK(NUM_node,TPP,NUM_max,TKK,TNN);

%%  determine whether the matrix coefficient is 0
        [TPP,TSS]=FMULTI(NUM_node,NNN,NN,TSS,TPP,TTT,TNN);

        for i = 1:NUM_upnode
            MRBPT0 = NODE_up(i,1);
            TKK(NN(MRBPT0)) = 1.0E35;
            TPP(MRBPT0) = NODE_up(i,2)*1.0E35;
        end
%%  GAUSSIAN  METHOD
        [TTT]=FGAUSSIAN(NUM_node,NNN,NN,TSS,TPP,TKK);
      
        T_NXNY=zeros(NY,NX);
        for i = 1:NY
            for j = 1:NX
                T_NXNY(i,j) =  TTT(NX*(NY-i)+j,1);
            end
        end
        T_T = T_NXNY;%
        clear T_NXNY
        AJG{day_run,1}=T_T; 
        clear c k
        str=['runing...',num2str(year_run),'year'];
        waitbar(day_run/365,h,str)
        day_total = day_total + 1;
    end
%%  save results
    nam = sprintf('%04d',year_run);
    na = [outpath,num2str(nam),'a.mat'];
    save(na,'AJG');
    clear AJG
end
close(h);
toc;



function [TTT]=FTTT(T_T,NX,NY)
%%      利用初始温度值建立温度场

%将建立的温度场分配序号
for i =  1:NY
    for j = 1:NX
        TTT(NX*(NY-i)+j,1) = T_T(i,j);
    end
end
end
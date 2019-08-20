%% 清空环境
clc;clear

%% 下载数据
%load scheduleData Jm T JmNumber
%工序 时间

%% 基本参数
NIND=40;        %个体数目
MAXGEN=500;      %最大遗传代数
GGAP=0.9;       %代沟
XOVR=0.8;       %交叉率
MUTR=0.2;       %变异率
gen=0;          %代计数器
%PNumber 工件个数 MNumber  工序个数
%[PNumber MNumber]=size(Jm);
PNumber=20;MNumber=2;JmNumber=8;
trace=zeros(2, MAXGEN);      %寻优结果的初始值
WNumber=PNumber*MNumber;     %工序总个数 80

%% 初始化
Number=zeros(1,PNumber);     % PNumber 工件个数
for i=1:PNumber
    Number(i)=MNumber;         %MNumber工序个数
end

% 代码2层，第一层工序，第二层机器
Chrom=zeros(NIND,2*WNumber);
for j=1:NIND
    WPNumberTemp=Number;
    for i=1:WNumber
        
        %随机产成工序
        val=unidrnd(PNumber);
        while WPNumberTemp(val)==0
            val=unidrnd(PNumber);
        end
        
        %第一层代码表示工序
        Chrom(j,i)= val;
        WPNumberTemp(val)=WPNumberTemp(val)-1;
        
        %第2层代码表示机器
        %随机产成工序机器
%        www = unidrnd(JmNumber);
        if mod(i,JmNumber)==0 %通过启发式思想生成顺序安排的机器，几乎一直是最优
            www=JmNumber;
        else
            www=mod(i,JmNumber);
        end
        Chrom(j,i+WNumber)=www;
        
    end
    for i=1:WNumber-1   %杜绝同一机器加工两次
        for k=i+1:WNumber
            if Chrom(j,i)==Chrom(j,k)
                while Chrom(j,i+WNumber)==Chrom(j,k+WNumber)
                    Chrom(j,k+WNumber)=unidrnd(JmNumber);
                end
            end
        end
    end
end
%计算目标函数值
[PVal ObjV P S]=cal(Chrom,PNumber,JmNumber); 
%% 循环寻找
while gen<MAXGEN
    
    %分配适应度值
    FitnV=ranking(ObjV);  %自动排序，给出最适合活下来的染色体
    %选择操作
    SelCh=select('rws', Chrom, FitnV, GGAP);       
    %交叉操作
    SelCh=across(SelCh,XOVR,PNumber,JmNumber);          
    %变异操作只变异机器
    SelCh=aberranceJm(SelCh,MUTR,PNumber,JmNumber);            
    
    %计算目标适应度值
    [PVal ObjVSel P S]=cal(SelCh,PNumber,JmNumber);   
    %重新插入新种群
    [Chrom ObjV] =reins(Chrom, SelCh,1, 1, ObjV, ObjVSel);       
    %代计数器增加
    gen=gen+1;       
    
    %保存最优值
    trace(1, gen)=min(ObjV);       
    trace(2, gen)=mean(ObjV);  
    
    % 记录最佳值
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);%最小时间
        STemp=S;
    end
    %记录 最小的工序
    if MinVal> trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
        STemp=S;
    end
    
end

% 当前最佳值
PVal=Val1; %工序时间
P=Val2;  %工序 
S=STemp; %调度基因含机器基因

%% 描绘解的变化
figure(1)
plot(trace(1,:));
hold on;
plot(trace(2,:),'-.');grid;
legend('解的变化','种群均值的变化');

%% 显示最优解
figure(2);
MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2);
for i=1:WNumber  
    val= P(1,i);
    a=(mod(val,100)); %工序
    b=((val-a)/100);  %工件
    Temp=[1:JmNumber];
    mText=Temp(MP(1,i));
    
    x1=PVal(1,i);
    x2=PVal(2,i);
    
    y1=mText-1;
    y2=mText;
    PlotRec(x1,x2,mText);
    
    PlotRec(PVal(1,i),PVal(2,i),mText);
    hold on;
    
    fill([x1,x2,x2,x1],[y1,y1,y2,y2],[1-1/b,1/b,b/PNumber]);
    text((x1+x2)/2,mText-0.25,num2str(P(i)),'FontSize',8);
end
    clc
    PVal
    MP
co=zeros(1,JmNumber);
usingtime=[16870        8980       17352        7792       16324        7792       15360        4870
];
for i=1:JmNumber
    for j=WNumber:-1:1
        if MP(j)==i
            co(i)=PVal(1,j); %给出结束时间索引
            break;
        end
    end
end
for i=1:JmNumber
    for j=WNumber:-1:1
        if MP(j)==i
            usingtime(i)=usingtime(i)+PVal(1,j)-PVal(2,j);
        end
    end
end
co
eo=MP(WNumber)
usingtime
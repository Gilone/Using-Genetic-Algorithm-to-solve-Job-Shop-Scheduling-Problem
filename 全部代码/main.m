%% ��ջ���
clc;clear

%% ��������
%load scheduleData Jm T JmNumber
%���� ʱ��

%% ��������
NIND=40;        %������Ŀ
MAXGEN=500;      %����Ŵ�����
GGAP=0.9;       %����
XOVR=0.8;       %������
MUTR=0.2;       %������
gen=0;          %��������
%PNumber �������� MNumber  �������
%[PNumber MNumber]=size(Jm);
PNumber=20;MNumber=2;JmNumber=8;
trace=zeros(2, MAXGEN);      %Ѱ�Ž���ĳ�ʼֵ
WNumber=PNumber*MNumber;     %�����ܸ��� 80

%% ��ʼ��
Number=zeros(1,PNumber);     % PNumber ��������
for i=1:PNumber
    Number(i)=MNumber;         %MNumber�������
end

% ����2�㣬��һ�㹤�򣬵ڶ������
Chrom=zeros(NIND,2*WNumber);
for j=1:NIND
    WPNumberTemp=Number;
    for i=1:WNumber
        
        %������ɹ���
        val=unidrnd(PNumber);
        while WPNumberTemp(val)==0
            val=unidrnd(PNumber);
        end
        
        %��һ������ʾ����
        Chrom(j,i)= val;
        WPNumberTemp(val)=WPNumberTemp(val)-1;
        
        %��2������ʾ����
        %������ɹ������
%        www = unidrnd(JmNumber);
        if mod(i,JmNumber)==0 %ͨ������ʽ˼������˳���ŵĻ���������һֱ������
            www=JmNumber;
        else
            www=mod(i,JmNumber);
        end
        Chrom(j,i+WNumber)=www;
        
    end
    for i=1:WNumber-1   %�ž�ͬһ�����ӹ�����
        for k=i+1:WNumber
            if Chrom(j,i)==Chrom(j,k)
                while Chrom(j,i+WNumber)==Chrom(j,k+WNumber)
                    Chrom(j,k+WNumber)=unidrnd(JmNumber);
                end
            end
        end
    end
end
%����Ŀ�꺯��ֵ
[PVal ObjV P S]=cal(Chrom,PNumber,JmNumber); 
%% ѭ��Ѱ��
while gen<MAXGEN
    
    %������Ӧ��ֵ
    FitnV=ranking(ObjV);  %�Զ����򣬸������ʺϻ�������Ⱦɫ��
    %ѡ�����
    SelCh=select('rws', Chrom, FitnV, GGAP);       
    %�������
    SelCh=across(SelCh,XOVR,PNumber,JmNumber);          
    %�������ֻ�������
    SelCh=aberranceJm(SelCh,MUTR,PNumber,JmNumber);            
    
    %����Ŀ����Ӧ��ֵ
    [PVal ObjVSel P S]=cal(SelCh,PNumber,JmNumber);   
    %���²�������Ⱥ
    [Chrom ObjV] =reins(Chrom, SelCh,1, 1, ObjV, ObjVSel);       
    %������������
    gen=gen+1;       
    
    %��������ֵ
    trace(1, gen)=min(ObjV);       
    trace(2, gen)=mean(ObjV);  
    
    % ��¼���ֵ
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);%��Сʱ��
        STemp=S;
    end
    %��¼ ��С�Ĺ���
    if MinVal> trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
        STemp=S;
    end
    
end

% ��ǰ���ֵ
PVal=Val1; %����ʱ��
P=Val2;  %���� 
S=STemp; %���Ȼ��򺬻�������

%% ����ı仯
figure(1)
plot(trace(1,:));
hold on;
plot(trace(2,:),'-.');grid;
legend('��ı仯','��Ⱥ��ֵ�ı仯');

%% ��ʾ���Ž�
figure(2);
MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2);
for i=1:WNumber  
    val= P(1,i);
    a=(mod(val,100)); %����
    b=((val-a)/100);  %����
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
            co(i)=PVal(1,j); %��������ʱ������
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
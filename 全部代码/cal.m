function [PVal ObjV P S]=cal(Chrom,PNumber,JmNumber)

% ����˵����       ���ݻ���Ⱥ,�������Ⱥ��ÿ������ĵ��ȹ���ʱ�䣬
%                 ������Сʱ��ĵ��ȹ���͵��ȹ���ʱ��
% ���������
%       Chrom     Ϊ������Ⱥ  
%       PNumber   ��������
%       JmNumber  ����number

% �������:
%       PVal      Ϊ��ѵ��ȹ���ʱ�� 
%       P         Ϊ�������ĵ��ȹ��� 
%       ObjV      ΪȺ��ÿ������ĵ��ȹ���ʱ��
%       S         Ϊ�������ĵ��Ȼ���

%��ʼ��
NIND=size(Chrom,1);
ObjV=zeros(NIND,1);

% ������� 
MNumber=2;

for i=1:NIND  
    
    %ȡһ������
    S=Chrom(i,:);
    
    %���ݻ��򣬼�����ȹ���
    P= calP(S,PNumber);
   
    %���ݵ��ȹ��򣬼�������ȹ���ʱ��
%    PVal=caltime(S,P,JmNumber,T,Jm); 
    PVal=operation2_1_1(S);    
    %ȡ���ʱ��
    MT=max(PVal);
    TVal=max(MT);  
    
    %����ʱ��
    ObjV(i,1)=TVal;
    
    %��ʼ��
    if i==1
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
    
    %��¼ ��С�ĵ��ȹ���ʱ�䡢��ѵ��ȹ���ʱ�� �������ĵ��ȹ���
    if MinVal> ObjV(i,1)
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end   
end 
 
%��ѵ��ȹ���ʱ�� �������ĵ��ȹ���
 PVal=Val1;
 P=Val2;
 S=STemp;

 
 
 
 
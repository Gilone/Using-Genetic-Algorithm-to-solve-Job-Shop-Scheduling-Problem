function ChromNew=aberranceJm(Chrom,MUTR,PNumber,JmNumber)

%��ʼ��
[NIND,WNumber]=size(Chrom);
WNumber=WNumber/2;

ChromNew=Chrom;
MNumber=2;
Number=zeros(1,PNumber);
for i=1:PNumber
  Number(i)=1;
end

for i=1:NIND    
                
    %ȡһ������
    S=Chrom(i,:);
                    
       for j=1:WNumber
          
            %�Ƿ����
          if MUTR>rand;
              
%               ѡ����������ѡ��
                S(j+WNumber)=unidrnd(JmNumber); 
          
                %ѡ������� �ӹ�ʱ���ٵ�ѡ���ʴ�
 %               if SizeTemp==1      
 %                      S(j+WNumber)=1;
 %               else
 %                   S(j+WNumber)=selectJm(S(j++WNumber),T{S(j),WPNumberTemp(S(j))});
 %               end
          end
          
        end         
   
  
    %���ݷ�����Ⱥ
    ChromNew(i,:)=S;
        
    for ii=1:WNumber-1   %�ž�ͬһ�����ӹ�����
        for kk=ii+1:WNumber
            if ChromNew(i,ii)==ChromNew(i,kk)
                while ChromNew(i,ii+WNumber)==ChromNew(i,kk+WNumber)
                    ChromNew(i,kk+WNumber)=unidrnd(JmNumber);
                end
            end
        end
    end
    
end

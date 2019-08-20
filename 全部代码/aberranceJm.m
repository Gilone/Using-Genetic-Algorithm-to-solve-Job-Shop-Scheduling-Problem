function ChromNew=aberranceJm(Chrom,MUTR,PNumber,JmNumber)

%初始化
[NIND,WNumber]=size(Chrom);
WNumber=WNumber/2;

ChromNew=Chrom;
MNumber=2;
Number=zeros(1,PNumber);
for i=1:PNumber
  Number(i)=1;
end

for i=1:NIND    
                
    %取一个个体
    S=Chrom(i,:);
                    
       for j=1:WNumber
          
            %是否变异
          if MUTR>rand;
              
%               选择机器（随机选择）
                S(j+WNumber)=unidrnd(JmNumber); 
          
                %选择机器（ 加工时间少的选择几率大）
 %               if SizeTemp==1      
 %                      S(j+WNumber)=1;
 %               else
 %                   S(j+WNumber)=selectJm(S(j++WNumber),T{S(j),WPNumberTemp(S(j))});
 %               end
          end
          
        end         
   
  
    %数据放入新群
    ChromNew(i,:)=S;
        
    for ii=1:WNumber-1   %杜绝同一机器加工两次
        for kk=ii+1:WNumber
            if ChromNew(i,ii)==ChromNew(i,kk)
                while ChromNew(i,ii+WNumber)==ChromNew(i,kk+WNumber)
                    ChromNew(i,kk+WNumber)=unidrnd(JmNumber);
                end
            end
        end
    end
    
end

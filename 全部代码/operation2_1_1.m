function PVal = operation2_1_1(K)

 M=8; %机器数
 N=40;%工序数
 n=K(1:40);
 m=K(41:80);
 T=[482 487 482 487 482 487 482 487;
    209 214 209 214 209 214 209 214];
%T=[428 431 428 431 428 431 428 431
%   406 409 406 409 406 409 406 409];
%T=[310 315 310 315 310 315 310 315
%    530 535 530 535 530 535 530 535];

 E=[27	27	45	45	59	59	73	73
32	32	50	50	64	64	78	78
45	45	27	27	45	45	59	59
50	50	32	32	50	50	64	64
59	59	45	45	27	27	45	45
64	64	50	50	32	32	50	50
73	73	59	59	45	45	27	27
78	78	64	64	50	50	32	32
]; %不清洗时间
 CO=[99999999       21325       24650       24195       24513       24695       24023       21588   
];
 EO=6; 
 C=zeros(1,N);
 PVal=zeros(2,N);
 fmindex=zeros(1,M);
 for j=1:M
    for i=1:N
        if m(i)==j
            fmindex(j)=i;
            break;
        end
    end
 end

i=1;
temp=[CO(m(i)) E(EO,m(i))];
C(i)=max(temp)+T(1,m(i));
PVal(2,i)=max(temp);
for i=2:N
    for k=i-1:-1:1 %判断当前加工的是第几步
        if n(k)==n(i)
            nstep=1;
            nmachine=k; %保留第一步工序的索引
        else
            nstep=0;
            nmachine=i;
        end
    end
    
    if isempty(find(fmindex==i, 1)) %判断机器之前有没有用过，先是已经用过的
        for j=i-1:-1:1
            if m(j)==m(i)  %找到上一次用这个机器，判断这个机器是否用完了
                x=j;
            break;                                           
            end
        end
        if nstep==0 % first step
            temp=[C(x) C(i-1)-T(m(i-1))+E(m(i-1),m(i))];
            C(i)=max(temp)+T(1,m(i));
            PVal(2,i)=max(temp);
        else
            if m(nmachine)==m(i-1)    %在此考虑三种绕路取物品的问题
                if nmachine~=i-1
                    temp=[C(x) C(i-1)-T(m(i-1))+E(m(i-1),m(i))]; %先被安排了另外的任务
                else
                    temp=[C(x) C(i-1)+E(m(i-1),m(i))]; %未被安排另外的任务
                end
            else
                temp=[C(x) C(i-1)-T(m(i-1))+E(m(i-1),m(nmachine))+E(m(nmachine),m(i)) C(nmachine)+E(m(nmachine),m(i))];
            end
            C(i)=max(temp)+T(2,m(i));
            PVal(2,i)=max(temp);
        end
    else  %这是第一次使用该机器
        
         if nstep==0 % first step
            temp=[CO(m(i)) C(i-1)-T(m(i-1))+E(m(i-1),m(i))];
            C(i)=max(temp)+T(1,m(i));
            PVal(2,i)=max(temp);            
        else
            if m(nmachine)==m(i-1)    %在此考虑三种绕路取物品的问题
                if nmachine~=i-1
                    temp=[CO(m(i)) C(i-1)-T(m(i-1))+E(m(i-1),m(i))]; %先被安排了另外的任务
                else
                    temp=[CO(m(i)) C(i-1)+E(m(i-1),m(i))]; %未被安排另外的任务
                end
            else
                temp=[CO(m(i)) C(i-1)-T(m(i-1))+E(m(i-1),m(nmachine))+E(m(nmachine),m(i)) C(nmachine)+E(m(nmachine),m(i))];
            end
            C(i)=max(temp)+T(2,m(i));
            PVal(2,i)=max(temp);     
         end
 
    end
end
PVal(1,:)=C;
end

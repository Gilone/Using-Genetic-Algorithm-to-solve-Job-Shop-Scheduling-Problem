function PVal = operation1_1_1(K)
%INPUT K IS AN ARRAY WHICH HAS 50 ELEMNTS,WHICH WE WANT
%M THE AMOUNT OF MACHINES N THE AMOUNT OF THINGS T THE ARRAY OF OPTIME 
%E THE MATRIAX OF WAYTIME
%TO THE ORIGINNAL NUMBER OF EACH MACHINE
 M=8;
 N=50;
 K=K(51:100);
%T=[572 577 572 577 572 577 572 577];
T=[588 591 588 591 588 591 588 591];
%T=[610 615 610 615 610 615 610 615];
 E=[53	53	73	73	86	86	99	99
56	56	76	76	89	89	102	102
73	73	53	53	73	73	86	86
76	76	56	56	76	76	89	89
86	86	73	73	53	53	73	73
89	89	76	76	56	56	76	76
99	99	86	86	73	73	53	53
102	102	89	89	76	76	56	56
];
 CO=[27146       28331       27275       27331       27404       27460       27533       28229
];
 EO=2;
 PVal=zeros(2,N);
 C=zeros(1,N);
 fmindex=zeros(1,M);
 for j=1:M
    for i=1:N
        if K(i)==j
            fmindex(j)=i;
            break;
        end
    end
 end

i=1;
temp=[CO(K(i)) E(EO,K(i))];
C(i)=max(temp)+T(K(i));
PVal(2,i)=max(temp);
for i=2:N
    if isempty(find(fmindex==i, 1))
        for j=i-1:-1:1
            if K(j)==K(i)
                x=j;
            break;
            end
        end
        temp=[C(x) C(i-1)-T(K(i-1))+E(K(i-1),K(i))];
        C(i)=max(temp)+T(K(i));
        PVal(2,i)=max(temp);
    else
        temp=[CO(K(i)),C(i-1)-T(K(i-1))+E(K(i-1),K(i))];
        C(i)=max(temp)+T(K(i));
        PVal(2,i)=max(temp);
    end
end
PVal(1,:)=C;

end


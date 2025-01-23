
a1vec=linspace(100,700,50);
a2vec=linspace(.1 ,20,50);
a3vec=linspace(.1,20,50);

    
            
xvec1=1:1:max_;
count=1;
tempval=Inf;

for i11=1:length(a1vec)
    a1param=a1vec(i11);
    for i22=1:length(a2vec)
    
        a2param=a2vec(i22);

        for i33=1:length(a3vec)
        
            a3param=a3vec(i33);
            temp=a1param*(xvec1.^(a2param-1)).*(exp(-xvec1/a3param))/(a3param^a2param)/gamma(a2param);
            tempval2= norm(thetahat(1:1:max_)-temp);
            
            if tempval2<tempval
            
                tempval=tempval2;
                a111=a1param; a222=a2param; a333=a3param;
                
            end
             
        end
        
    end
    
end
            
thetahatw=a111*(xvec1.^(a222-1)).*(exp(-xvec1/a333))/(a333^a222)/gamma(a222);
[a111 a222 a333]


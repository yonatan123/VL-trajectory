% Generate samples using the multivariate normal model. 
%INPUT: n- number of samples
%       yt- mean ct-value (d-vector)
%       p0- prior probability of day of 1st measurement
%       Q-  covariance (d x d positive definite matrix)
%OUTPUT: structure D with fields: ct2: a matrix with each row consisting a
%        pair of ct-value measurements, and delta2: a column vector with
%        the corresponding delta's
function D=generatesamples(n,yt,p0,Q) 

    max_=length(yt);
    
    x1=[];, delta=[];, y1=[];, y2=[]; 
    
    while 1
        %sample integer from p0
        cp = [0, cumsum(p0)];
        r = rand;
        temp1 = find(r>cp, 1, 'last'); 
        
        % samples Delta between two ct measurements from Unif(2,3)
        dtemp=1+randi(2);
        temp2=temp1+dtemp;
        
        x1=[x1 temp1]; delta=[delta dtemp];  
        
        tempcov=Q([temp1  temp2],[temp1  temp2]);
       
        %sample ct-value from the two sampled days using multivariate
        %normal with parameters yt, Q
        tempz=[yt(temp1) yt(temp2) ]'+sqrtm(tempcov)*randn(2,1);
        y1=[y1 tempz(1)]; y2=[y2 tempz(2)];
        
     
    
       if length(x1)==n
           break
        end
        
    end
    
    
    
    [a b]=sort(delta);
    delta=delta(b);, y1=y1(b); y2=y2(b); x1=x1(b);


    delta1=delta; ind=find(delta>14); y1(ind)=[]; y2(ind)=[]; x1(ind)=[]; delta(ind)=[];
    %gather all samples in a struct D
    D=struct; D.ct2=[y1' y2']; D.delta2=delta';

    figure(100), scatter(x1,y1,'.b'); hold on;, scatter(x1+delta,y2,'.r'); %plot(yt,'k','LineWidth',1.5);

end
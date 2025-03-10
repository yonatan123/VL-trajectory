 


%EM-algorithm: estimation under unimodality constraints
%Input: D: data
%       theta0, sigma0, q0: intial mean, covariance and q
%       maxiter, tol: maximum number of iterations, tolerance
%       maxpoint: index of unimodal peak
%Output:thetahat, sigmahat, qhat: estimated mean, covariance and q
%       lvec: vector of marginal likelihood values, computed for each iteration
function [thetahat sigmahat qhat   lvec]=emfitunimodal(D,theta0,sigma0,q0,maxiter,tol,maxpoint)
        
        %figure, hold on
        options = optimoptions('quadprog','Display','off');
        dimq=size(sigma0,1);
        max_=length(q0);
        thetahat=theta0;
        qhat=q0;
        sigmahat=sigma0;

        error=1;
    
        lvec=[Inf];
        
        
        if maxpoint>0
            Aineq=eye(dimq)-diag(ones(dimq-1,1),1);
            Aineq=Aineq(1:dimq-1,:);
            Aineq(maxpoint:end,:)=-Aineq(maxpoint:end,:);   
        end
        Aeq1=zeros(dimq-max_,dimq);
        Aeq1(:,max_)=1;
        Aeq1( 1:end ,max_+1:end)=-eye(dimq-max_);
        

        numgroups=length(fieldnames(D))/2;
        
        
        for loop=1:maxiter
          
            if mod(loop,50)==1, 
                strtemp=strcat(['Iteration:',num2str(loop) ]);, 
                disp(strtemp)
            end            

            thetaprev=thetahat;
            sigmaprev=sigmahat;
            qprev=qhat;
                       


            thetahatcurrent=0*thetahat; 
            thetahatcurrentnoe=0*thetahat; 


       


            
            count=0;
            numobs=1;
            ntotal=0;
            
            numerator=0;
            denominator=0;

            val=0;

            while count<numgroups
                

                numobs=numobs+1;
                try
                    y=D.(strcat(['ct',num2str(numobs)]));
                    delta=D.(strcat(['delta',num2str(numobs)]));
                catch
                    continue
                end
                count=count+1;


                dim=size(y,2);
                n=size(y,1);
                ntotal=ntotal+n;
                cdim=dimq-dim;
                %etemp=[];
                %Etheta=[];
               
                
                  
                %calculate E
                E1=zeros(n,max_);
                for i=1:max_
                    
                    for j=1:n
    
                        temp=[i i+delta(j,:)];
                        S=sigmahat(temp,temp);
                        M=thetahat(temp);    
                        E1(j,i)=mvnpdf(y(j,:),M,S)*qhat(i);
                    end
    
                end    


                
                val=val+sum( log(sum(E1,2))  );
                E1=E1./repmat(sum(E1,2),1,max_);
                
                numerator=numerator+sum(E1,1);
                denominator=denominator+sum(sum(E1,1));




                %calculate Ez and Ezz'
                zlatent1=zeros(dimq,1);  
                covlatent1=zeros(dimq);
            
            
                  for j=1:n
                    
                      for i=1:max_
                        
                          index=[i i+delta(j,:)];
                          cindex=1:1:dimq;
                          cindex(index)=[];
                          temp=thetahat(cindex)'+  sigmahat(cindex,index)*inv(sigmahat(index,index))*(y(j,:)'-thetahat(index)');
                          eij=zeros(dimq,1);
                          eij(cindex)=temp; eij(index)=y(j,:)';
                          
                          
                          %zlatent
                          zlatent1(cindex)=zlatent1(cindex)+temp*E1(j,i);
                          zlatent1(index)=zlatent1(index)+y(j,:)'*E1(j,i);
            
            
                          %covlatent
                          temp=sigmahat(cindex,cindex) - sigmahat(cindex,index)*inv(sigmahat(index,index))*sigmahat(index,cindex);
                          C=eij*eij';
                          C(cindex,cindex)=C(cindex,cindex)+temp;
                          C=C - thetahat'* eij' - (thetahat'* eij')'+ (thetahat'* thetahat);
                          covlatent1 = covlatent1 + E1(j,i)*C; 
            
            
            
                      end
            
                  end








            end
            lvec=[lvec val];
            

            %%%%%%%%%%%%%%%%%%%%%%%%


            thetahat=zlatent1'/ntotal;
            %sigmahat=covlatent/ntotal;
            qhat=numerator/denominator;
            

            H=1*inv(sigmahat); f=inv(sigmahat)*thetahat'; 
            if maxpoint>0
                [thetahat tempval]=quadprog(H,-f,Aineq,zeros(dimq-1,1),Aeq1,zeros(size(Aeq1,1),1),zeros(dimq,1),[],[],options);, 
            end
            

            thetahat=thetahat';
            %plot(thetahat); pause(.1)


            p=100;
            rhovar=linspace(0,.99,p);
            %rhovar=[.1]; p=length(rhovar);
            objectivevalue=Inf;
            for j=1:p
                rho=rhovar(j);
                sigmatemp=toeplitz([rho.^(0:1:dimq-1)]);
                tracetemp=trace(inv(sigmatemp)*covlatent1) ;
                s= tracetemp/dimq/ntotal;
                valuetemp=dimq*ntotal*log(s)+ntotal*log(det(sigmatemp))+tracetemp/s;
                %valuetemp=ntotal*log(det(s*sigmatemp))+tracetemp/s;
                if valuetemp<objectivevalue
                    current=[s rho];
                    objectivevalue=valuetemp;
                end
            end
            sigmahat=current(1)*toeplitz([current(2).^(0:1:dimq-1)]);
            

            error=norm([sigmahat(1,1) sigmahat(1,2)] -[sigmaprev(1,1) sigmaprev(1,2)])^2+norm(qhat-qprev)^2+norm(thetahat-thetaprev)^2;
            
            if sqrt(error)<tol
                break
            end

        end


end



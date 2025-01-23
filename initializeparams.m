

        theta0=rand(1,max_);
        theta0(max_:2*max_)=theta0(end);
                        
        inits=2*rand(2,1);
        s=max(inits); rho=min(inits);
        sigma0=toeplitz([s rho.^(1:1:2*max_-1)]);
                            
        q0=rand(1,max_);
        q0=q0/sum(q0);
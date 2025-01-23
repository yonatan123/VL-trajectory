
%run synthetic trial. Estimation using synthetic samples, assuming no
%constraints, unimodal and Gamma
clear, close all, clc

load('synthparameters.mat'); %load true parameters
Q=Q/100;  %scale covariance by 1/100
n=500;    % #samples
D=generatesamples(n,yt,p0,Q);
pause(.1)

maxiter=1000;
tol=.01;
numtrials=5;
countrials=0;



%estimate without constraints, with 'numtrials' number of random intializations
disp('Estimation without constraints')
disp('---------------------------------------------------------')
likelihoodvalue=-inf;
likelihoodvector=[];

while countrials<=numtrials
    try
        initializeparams
        
        [thetahat sigmahat qhat    lvec]=emfit(D,theta0,sigma0,q0,maxiter,tol);
        
        if lvec(end)>likelihoodvalue
            likelihoodvector=[likelihoodvector lvec(end)];
            likelihoodvalue=lvec(end);
            thetahatfinal=thetahat;, sigmahatfinal=sigmahat;, qhatfinal=qhat;    
        end

        countrials=countrials+1;
        disp(strcat(['# trial: ', num2str(countrials)]))
        


    catch disp('initialization error')
    end
end
figure(100),  plot(thetahatfinal, 'g','LineWidth',1.5)            
pause(.1)



disp('Estimation under unimodality constraint')
disp('---------------------------------------------------------')
%unimodal estimate 
likelihoodvalue=-inf;
likelihoodvector=[];
countrials=0;

while countrials<=numtrials
    maxpoint=2;
    while maxpoint<7
        try
            initializeparams
            
            [thetahat sigmahat qhat    lvec]=emfitunimodal(D,theta0,sigma0,q0,maxiter,tol,maxpoint);
            
            if lvec(end)>likelihoodvalue
                likelihoodvector=[likelihoodvector lvec(end)];
                likelihoodvalue=lvec(end);
                thetahatfinal=thetahat;, sigmahatfinal=sigmahat;, qhatfinal=qhat;    
            end


            disp(strcat(['# trial: ', num2str(countrials), ' maxpoint: ' , num2str(maxpoint)]))    
            maxpoint=maxpoint+1;


    
        catch disp('initialization error')
        end
        

    end
    countrials=countrials+1;
end
figure(100),  plot(thetahatfinal, 'm','LineWidth',1.5)            
pause(.1)





disp('Estimation under Gamma constraint')
disp('---------------------------------------------------------')
%Gamma estimate 
likelihoodvalue=-inf;
likelihoodvector=[];
countrials=0;

while countrials<=numtrials
    maxpoint=2;
    while maxpoint<7
        try
            initializeparams
            
            [thetahat sigmahat qhat    lvec]=emfitunimodal(D,theta0,sigma0,q0,maxiter,tol,maxpoint);
                        
            gamma_approximation
            
            a1vec=linspace(.75*a111,1.25*a111,50);
            a2vec=linspace(.75*a222,1.25*a222,50);
            a3vec=linspace(.75*a333,1.25*a333,50);
            [thetahatg sigmahatg qhatg   lvec]=emfitgamma(D,thetahat,sigmahat,qhat,maxiter,tol,a1vec,a2vec,a3vec);
                                
            if lvec(end)>likelihoodvalue
                likelihoodvector=[likelihoodvector lvec(end)];
                likelihoodvalue=lvec(end);
                thetahatfinal=thetahatg;, sigmahatfinal=sigmahatg;, qhatfinal=qhatg;    
            end


            disp(strcat(['# trial: ', num2str(countrials), ' maxpoint: ' , num2str(maxpoint)]))    
            maxpoint=maxpoint+1;


    
        catch disp('initialization error')
        end
        

    end
    countrials=countrials+1;
end



figure(100),  plot(thetahatfinal, 'c','LineWidth',1.5)  
figure(100), legend('Ct1 samples','Ct2 samples','no constraints','unimodal','Gamma')

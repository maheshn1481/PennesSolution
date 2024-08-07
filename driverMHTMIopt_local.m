clear all
clc

myoptions.Algorithm = 'constDirect'
%myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-9, 'OptimalityTolerance',2.5e-4,'Algorithm','interior-point','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'HonorBounds',true, 'Diagnostic','on','FunValCheck','on' )

%% % monitor memory: while [ -e /proc/3291925 ] ; do  top -b -n 1 -p 3291925 >>process.txt ;sleep 60; done

%function driverHPMIopt(NGauss,NumberUncertain,modelSNR,myoptions,ObjectiveType,GaussLegendre )
NGauss = 3,NumberUncertain=1,modelSNR=10, ObjectiveType = 'TotalSignal',GaussLegendre=false


NGauss,NumberUncertain,modelSNR,myoptions.Algorithm,ObjectiveType,GaussLegendre
close all
%% Tissue Parameters
kmean = [ .5 ]; %  Thermal Conductivity
kstdd = [ .2 ]; % s
klb   = [ 5  ]; % s
kub   = [ 45 ]; % s
wmean = [ 6  ]; % Blood perfusion rate
wstdd = [ 3  ]; % s
wlb   = [ 2  ]; % s
wub   = [ 10 ]; % s
crhomean = [ 30000 ];       % Density, Specific Heat Production
crhostdd = [ 10000 ];       % s
crholb   = [ 20000 ];       % s
crhoub   = [ 40000 ];       % s
chimean = [ 0.05 ];      % What is this?
chistdd = [ .01  ];      % s
chilb   = [ 0.01 ];      % s
chiub   = [ 0.20 ];      % s
ubmean  = [ 4    ];      % What is this?
ubsttd  = [ 1.3  ];      % s
ublb    = [ 0    ];       % s
ubub    = [ 7    ];       % s

%% Variable Setup
Ntime = 30;
currentTR = 3;
timelist = (0:(Ntime-1))*currentTR ;
Power = 20*ones(Ntime,1);
temperature = zeros(Ntime,1);

for iii = 1:Ntime
    temperature(iii) = pennesmht(20,kmean,wmean,Power(iii),.001);
end

%% Plot initial guess
plotinit = true;
if plotinit
    % plot initial guess
    figure(1)
    plot(timelist ,temperature)
    ylabel('temperature ')
    xlabel('sec')
end


%% optimize MI
optf = true;
if optf

    tic;
    % setup optimization variables
    Nspecies = 1
    powerList = optimvar('powerList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35);

    signu = sqrt(2* Ntime) * 10;
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));

    % switch between uniform and Gaussian RV for prior
    if(GaussLegendre)
        QuadratureRule = 'Legendre';
    else
        QuadratureRule = 'Hermite';
    end
    switch (NumberUncertain)
        case(1)
            if(GaussLegendre)
                [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,chilb,chiub);
            else
                [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,chimean,chistdd);
            end
            kqp     = kmean;
            wqp     = wmean;
            crhoqp  = crhomean;
            chiqp   = xn{1}(:);
        case(4)
            if(GaussLegendre)
                [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,[klb; wlb; crholb; chilb],[kub; wub; crhoub; chiub]);
            else
                [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[kmean; wmean; crhomean; chimean],[kstdd; wstdd; crhostdd; chistdd]);
            end
            kqp    = xn{1}(:);
            wqp    = xn{2}(:);
            crhoqp = xn{3}(:);
            chiqp  = xn{4}(:);
    end
    %alphaqp = xn{6}(:);
    %betaqp  = xn{7}(:);

    lqp=length(xn{1}(:));
    %statevariable    = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0,'UpperBound',.1);
    statevariable  = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0);
    auxvariable      = optimexpr(   [Ntime,Nspecies,lqp]);
    stateconstraint  = optimconstr(    [Ntime,Nspecies,lqp]);

    disp('build state variable')

    % [expATRoneone(end),0; expATRtwoone(end), expATRtwotwo(end)]
    % IC
    stateconstraint(1,:,:)  = statevariable(1,:,:) ==0;
    auxvariable(1,:,:) =0;
    
    for iii = 1:Ntime-1
        % setup state as linear constraint
        auxvariable(iii+1,1,:) =  pennesmht(kqp,wqp,crhoqp,chiqp,powerList(iii+1));
        stateconstraint(iii+1,1,:)  = statevariable(iii+1,1,:) ==  pennesmht(kqp,wqp,crhoqp,chiqp,powerList(iii+1));
    end


    disp('build objective function')
    expandvar  = ones(1,lqp);

    switch (ObjectiveType)
        case('TotalSignal')
            % TODO - repmat does not work well with AD
            % TODO - replace repmat with matrix
            % sumstatevariable = squeeze(sum(repmat(sin(FaList)',1,1,lqp).*statevariable,1));
            sumstatevariable = optimexpr([Nspecies,lqp]);
            for jjj = 1:lqp
                sumstatevariable(:,jjj) =  sum(statevariable(:,:,jjj) ,1)';

            end
            diffsumm =sumstatevariable(1,:)' * expandvar   - expandvar' * sumstatevariable(1,:);
            negHz = 0;
            for jjj=1:lqp2
                znu=xn2{1}(jjj) ;
                % note that the sqrt(pi)^N+1  and 2^N factors from the integration over priors is included in the quadrature routines.
                % this makes is easier to switch between Gaussian and Uniform RV
                negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
            end
    end
    % MI = H(z) - H(z|P)
    %  H(z|P)  constant ==> max MI = max H(z) = min -H(z)
    MIGaussObj = negHz;

    %%
    % Create an optimization problem using these converted optimization expressions.
    disp('create optim prob')
    convprob = optimproblem('Objective',MIGaussObj , "Constraints",stateconstraint);
    myidx = varindex(convprob )
    %%
    % View the new problem.
    %show(convprob)
    problem = prob2struct(convprob,'ObjectiveFunctionName','reducedObjective','ConstraintFunctionName','reducedConstraint');


    % compare walker solution at qp
    switch (NumberUncertain)
        case(3)
            optparams = struct('t0',[t0qp(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
                'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplqp(1) ;0,0],...
                'TRList',TR_list,'PerfusionTerms',[kveqp(1),0],'volumeFractions',ve,...
                'fitOptions', opts);
        case(4)
            error("WIP")
        case(5)
            error("WIP")
    end
    % solve
    switch (myoptions.Algorithm)
        case('constDirect')
            powergrid = [0:1:35];
            brutesearch= zeros(size(powergrid ));
            backspaces = '';
            for iii = 1:length(powergrid(:))
                %disp(sprintf('%d %d ',pyrgrid(iii),lacgrid(iii)));
                % Print percentage progress
                percentage = iii/length(powergrid(:));
                perc_str = sprintf('completed %3.1f', percentage);
                fprintf([backspaces, perc_str]);
                backspaces = repmat(sprintf('\b'), 1, length(perc_str));

                x0.powerList = repmat(powergrid(iii),1,Ntime);
                x0.state  = evaluate(auxvariable ,x0);
                brutesearch(iii) = evaluate(MIGaussObj,x0);
            end
            save(sprintf('brutesearchNG%dNu%d%s%sSNR%02d%s.mat',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule) ,'brutesearch','powergrid')
            [maxMI,idmax] = max(brutesearch(:));
            [fval,idmin] = min(brutesearch(:));
            popt.powerList = repmat(powergrid(idmin),1,Ntime);
            %pneg.FaList = repmat([pi/180*pyrgrid(idmax);pi/180*lacgrid(idmax)],1,Ntime);
            %optparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);
            %refparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);

            handle = figure(5)
            plot(powergrid,brutesearch)
            xlabel('power ')
            ylabel('MI ')
            title(sprintf('max %f min %f',maxMI,fval) )
            %text(pyrgrid(idmin)+1,lacgrid(idmin)+1, sprintf('opt %d %d', pyrgrid(idmin), lacgrid(idmin)));
            %text(pyrgrid(idmax)+1,lacgrid(idmax)+1, 'control');
        otherwise
            InitialGuess =  [flips(:)];
            pmin =  [flips(:)*0];
            pmax =  [flips(:)*0+35*pi/180];
            tolx=1.e-9;
            tolfun=5.e-4;
            maxiter=400;

            % truthconstraint = infeasibility(stateconstraint,x0);
            %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ConstraintDerivative', 'auto-reverse', 'ObjectiveDerivative', 'auto-reverse' )
            Fx = @(x) MIGHQuadHPTofts(x, problem, myidx,Nspecies,Ntime,auxvariable);
            [designopt,fval,exitflag,output,lambda,grad,hessian] ...
                =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],myoptions);

            handle = figure(5)
            optparams.FaList = reshape(designopt(:),size(params.FaList ));
            refparams.FaList = reshape(designopt(:),size(params.FaList ));
            popt.FaList      = reshape(designopt(:),size(params.FaList ));
    end
    % save convergence history
    set(gca,'FontSize',16)
    saveas(handle,sprintf('historyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule ),'png')
    popt.state       = evaluate(auxvariable, popt);
    toc;


    %% [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',optparams);
    %% [t_axisref,Mxyref,Mzref] = model.compile(M0.',refparams);
    %% save(sprintf('poptNG%dNu%d%s%sSNR%02d%s.mat',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule) ,'fval','popt','params','Mxy','Mz','Mxyref','Mzref','signu','signuImage')
    %% handle = figure(10)
    %% plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
    %% ylabel('MI Mxy')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptMxyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
    %% handle = figure(11)
    %% plot(params.TRList,popt.FaList(1,:)*180/pi,'b',params.TRList,popt.FaList(2,:)*180/pi,'k')
    %% ylabel('MI FA (deg)')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% ylim([0 40])
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptFANG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
    %% handle = figure(12)
    %% plot(params.TRList,Mzopt(1,:),'b--',params.TRList,Mzopt(2,:),'k--')
    %% hold
    %% plot(params.TRList,cos(optparams.FaList(1,:))'.*popt.state(:,1, 1),'b',params.TRList,cos(optparams.FaList(2,:))'.*popt.state(:,2, 1),'k')
    %% %if(lqp > 1)
    %% %  plot(params.TRList,popt.state(:,1, 5),'b',params.TRList,popt.state(:,2, 5),'k')
    %% %  plot(params.TRList,popt.state(:,1,10),'b',params.TRList,popt.state(:,2,10),'k')
    %% %  plot(params.TRList,popt.state(:,1,15),'b',params.TRList,popt.state(:,2,15),'k')
    %% %end
    %% ylabel('MI Mz ')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptMzNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
end


function [MIobjfun, MIobjfun_Der]=MIGHQuadHPTofts(xopt,problem,myidx,Nspecies,Ntime,auxvariable)
x0.FaList = reshape(xopt,Nspecies,Ntime);
x0.state  = evaluate(auxvariable ,x0);
Xfull = [ x0.FaList(:); x0.state(:)];
[MIobjfun,initVals.g] = problem.objective(Xfull);
[initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
objectiveGradFA    = initVals.g(myidx.FaList);
objectiveGradState = initVals.g(myidx.state);
jacobianFA    = initConst.ceqGrad(myidx.FaList,:);
jacobianState = initConst.ceqGrad(myidx.state,:);
adjointvar =-jacobianState \objectiveGradState ;
MIobjfun_Der = objectiveGradFA +  jacobianFA *   adjointvar ;
end
function tempqoi = pennesmht(k,w,crho,chi,P)
mua   = 0.45e+2
mus   = 47.0e+2
anfact= .9
mutr  = mua+mus*(1.0-anfact)
mueff = sqrt(3.0*mua*mutr)
ua    =  310
R1    =  .001
R2    =  .03
cblood = 3840.0
r = .001
u0 = 21;
s1 = 3.0/4.0/pi*P*mua*mutr/(w-k*mueff*mueff)*exp(-mueff*r)/r+ua;      s2 = s1;      s5 = 1/r*exp(sqrt(w/k)*r)*(-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*pi*R1*w+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+3.0*sqrt(w/k)*R2*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*pi*R1*w-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff-3.0*P*mua*mutr*mueff*R2*exp(-mueff*R2-sqrt(w/k)*R1)-3.0*P*mua*mutr*exp(-mueff*R2-sqrt(w/k)*R1)+4.0*exp(-sqrt(w/k)*R2)*ua*pi*R1*w-4.0*exp(-sqrt(w/k)*R2)*u0*pi*R1*w+4.0*exp(-sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+3.0*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)-4.0*exp(-sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff)/4.0;
s6 = exp(-sqrt(w/k)*(-R1+R2))/(-w+k*mueff*mueff)/pi/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);      s4 = s5*s6;      s6 = 1/r*exp(-sqrt(w/k)*r)*exp(-sqrt(w/k)*(-R1+R2))/4.0;
s9 = 4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*w-4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*w*sqrt(w/k)*R2-4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff*sqrt(w/k)*R2-3.0*P*mua*mutr*exp(sqrt(w/k)*R2-mueff*R1)+3.0*P*mua*mutr*sqrt(w/k)*R2*exp(sqrt(w/k)*R2-mueff*R1)-4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*w+4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*w*sqrt(w/k)*R2+4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff-4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff*sqrt(w/k)*R2+3.0*P*mua*mutr*mueff*R2*exp(sqrt(w/k)*R1-mueff*R2)+3.0*P*mua*mutr*exp(sqrt(w/k)*R1-mueff*R2);      s10 = 1/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);      s8 = s9*s10;      s9 = 1/pi/(-w+k*mueff*mueff);      s7 = s8*s9;      s5 = s6*s7;      s3 = s4+s5;
tempqoi  = s2+s3;
end
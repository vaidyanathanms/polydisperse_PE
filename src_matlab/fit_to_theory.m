% To fit to function hexp(beta EN)/(1+h exp(beta EN))
function [hval,eval] = fit_to_theory(xdata,ydata,beta,c0)
fitfun = @(x,xdata)((x(1)*exp(x(2)*beta*xdata))./(1+x(1)*exp(x(2)*beta*xdata)));
lb = [0,-Inf]; ub = [1,Inf]; %lower/upper bounds
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
%[bad_guess] = lsqcurvefit(fitfun,c0,xdata,ydata,lb,ub,options);
[coeff,resnorm] = lsqcurvefit(fitfun,c0,xdata,ydata,lb,ub,options);
hval  = coeff(1); eval = coeff(2);
fprintf('hval: %g\t eval: %g\t resnorm: %g\n', hval,eval,resnorm)


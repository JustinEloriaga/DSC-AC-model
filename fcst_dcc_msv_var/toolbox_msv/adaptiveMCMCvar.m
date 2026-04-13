function logvarnew = adaptiveMCMCvar(logvarold, alp, alphat, b, c, simulind)
% adaptive mcmc - proposal log(variance)
% logvarold: previous proposal log variance
% alp: mhratio today
% alphat: mhratio target
% b: bound
% c: cooling rate
% simulind: iteration index

% note: when c=inf, no adaptation

logvarnew = real(logvarold + ((simulind+1)^(-c))*(alp-alphat));

if logvarnew < -b
    logvarnew = -b;
elseif logvarnew > b
    logvarnew = b;
end
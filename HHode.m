
function xdot= HHode(t, x, flag, transitions, const)  % Added transitions
% ODE variables
V=x(1); m=x(2); h=x(3); n=x(4);

% Calculates all the alphas and beta
nAlpha = 0.01*(V - (const.vRest + 10))/(1 - exp(-(V - (const.vRest + 10))/10));
nBeta = 0.125*exp(-(V - const.vRest)/80);
mAlpha = 0.1*(V - (const.vRest + 25))/(1 - exp(-(V - (const.vRest + 25))/10));
mBeta = 4*exp(-(V - const.vRest)/18);
hAlpha = 0.07*exp(-(V - const.vRest)/20);
hBeta = 1/(1 + exp(-(V - (const.vRest + 30))/10));

% Calculates tau and inf
tauN = 1/(nAlpha + nBeta);
infN = nAlpha * tauN;

tauM = 1/(mAlpha + mBeta);
infM = mAlpha * tauM;

tauH = 1/(hAlpha + hBeta);
infH = hAlpha * tauH;

% Input current function generator
if ~isempty(transitions)
    Iapp = @(t) 10 * exp(-((t - transitions(1))*2).^2); % Initialize Iapp
    for i = 2 : length(transitions)
        % Add more inputs
        Iapp = @(t) Iapp(t) + 10 * exp(-((t - transitions(i))*2).^2);
    end
else
    Iapp = @(t) 0*t; % Set input to 0 when there is no transition 
end

% calculates the ODE
xdot(1,1) = (Iapp(t) - const.gNa * m^3 * h * (V - const.eNa) - const.gK * (V - const.eK) * n^4 - const.gLeak * (V - const.eLeak))/const.C;
xdot(2,1) = -(m - infM)/tauM;
xdot(3,1) = -(h - infH)/tauH;
xdot(4,1) = -(n - infN)/tauN;
end
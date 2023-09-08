function para = atmos(H, ID)

    T0 = 288.15;
    p0 = 101325;
    rho0 = 1.225;
    
    %lapse rate constant
    kappa = -0.0065 / T0;

if (H < 11000)
        theta = 1 - 0.000022558 * H;
        TempRatio = theta;
        p = p0 * theta .^ 5.2561;
        PressRatio = p ./ p0;
        rho = rho0 * theta .^ 4.2561;
        DensRatio = rho ./ rho0;
        SoundSpeed = sqrt(1.4 * 286.9 * TempRatio * T0);
      
elseif ((11000<H)<65617/3.281)
    TempRatio = 0.751865;
    PressRatio = 0.223361*exp(-(H*3.281-36089)/20806);
    DensRatio = 0.297076*exp(-(H*3.281-36089)/20806);
elseif((65617/3.281<H)<104987/3.281)
    TempRatio = 0.682457 + H*3.281/945374;
    PressRatio = (0.988626+H*3.281/652600)^(-34.1632);
    DensRatio = (0.978261+H*3.281/659515)^(-35.16320);
elseif ((104987/3.281<H)<154199/3.281)
    TempRatio = 0.482561 + H*3.281/337634;
    PressRatio = (0.898309 + H*3.281/181373)^(-12.20114);
    DensRatio = (0.857003+H*3.281/190115)^(-13.20114);
elseif((154199/3.281<H)<167323/3.281)
    TempRatio = 0.939268;
    PressRatio = exp(0.00109456*exp(-(H*3.281-154199)/25992));
    DensRatio = 0.00116533*exp(-(H*3.281 - 154199)/25992);
elseif((167323/3.281<H)<232940/3.281)
    TempRatio = 1.434843 - H*3.281/337634;
    PressRatio = (0.838263-H*3.281/577922)^(12.20114);
    DensRatio = (0.79899-H*3.281/606330)^(11.20114);
elseif((232940/3.281<H)<278386/3.281)
    TempRatio = 1.237723 - H*3.281/472687;
    PressRatio = (0.917131-H*3.281/637919)^(17.0816);
    DensRatio = (0.900194 - H*3.281/649922)^(16.08160);
end
    
   switch ID
    case 0
        para = TempRatio;
    case 1
        para = PressRatio;
    case 2
        para = DensRatio;
    case 10
        para = TempRatio * T0;
    case 12
        para = DensRatio * rho0;
    case 13
        para = SoundSpeed;
   end

end

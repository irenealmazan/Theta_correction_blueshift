%scale U to be 0 <= U <= 1 
U = U / max(max(U)) + 1;
U = U / max(max(U));

%See "Magneto-optics of Gd and Tb in the soft x-ray resonance regions"
%paper, figure 3 for Gd

%In U(x,y), I'll call 1 --> ferromagnet, 0 --> antiferromagnet
%For ferromagnet, i.e.      U(x,y) = 1      mu [nm^-1] = 0.12
%For antiferromagnet,i.e.   U(x,y) = 0      mu [nm^-1] = 0.07

mu_ferromagnet = 0.12*1000; %in um^-1
mu_antiferromagnet = 0.07*1000; %in um^-1

%so, to calculate attenuation, need to create an array using U that
%instead of being scaled from 0 to 1 is scaled from 0.07 to 0.12, multiply
%by thickness of Gd, multiply by -1 and take exponential ?

mu_Gd = U*(mu_ferromagnet - mu_antiferromagnet) + mu_antiferromagnet;
d_Gd = 0.04545;                                 %in [um]
U_Gd_Attenuation = U .* exp(-1 * mu_Gd * d_Gd);

%Attenuation in Fe is not dependent on the type of domain (we're tuned to
%Gd edge at Ian's beamline) i.e. there's no dichroism between x-rays and Fe

mu_Fe = 0.15*1000;% in [um]; check to make sure this is correct at 1.189 keV 
d_Fe = 0.04141;                                 %in [um]
U_Fe_Attenuation = U * exp(-1 * mu_Fe * d_Fe);

U_Fe_Gd_Attenuation = U_Fe_Attenuation .* U_Gd_Attenuation;

%Nov 1 2008...looks like that with such a small intensity in the U_Fe_Gd_Attenuation 
%matrix, the algorithm is having trouble converging. I'll try to scale it
%by 1E5 or so, and see what happens...

U_Fe_Gd_Attenuation = U_Fe_Gd_Attenuation * 1E5;

function [out_mu, out_rho, out_fail] = absorb( e_name, energy, thick)% [mu, rho, ifail] = absorb( name, energy, thick)% Calculate absorption coefficient mu (1/microns), density (gms/cm^3)% from e_name(string, either element, material or crystal)% energy (eV), and thick ( microns).%   calls to: anomal, atomdata, element, ff_to_mu, gtline, matlabs, %   rdlat, struct_f% This can be called in two ways, either as a function with returned% variables or as a script(-like) without output variables, at which % point certain info is printed on the screen.  To wit:% >> [mu,rho]=absorb(42,3000,3)% mu = 2.0321% rho = 10.2000% or% >> absorb si 8048% Element Si, Z= 14  Amu=     28.1  Rho=     2.32 gms/cc% f'=  0.254 f"=   0.33% Rayleigh= 0.00424 Compton=  0.00052% Absorption length=   69.1 microns at   8048 eV% Transmitted Intensity=  0.235 for    100 microns% Fortran version written by Sean Brennan at SSRL, described in% S. Brennan and P.L. Cowan, Rev. Sci. Instrum., 63, 850 (1992).% ported to matlab by Alex Lessman, SMB & Anneli Munkholmif nargin == 0,	e_name= gtline('Element name','Si','s');endif nargin < 2,	energy= gtline('Energy [eV]',8000);else	if isstr(energy),		energy= str2num(energy);	endendif nargin < 3,	thick= 100;		% defaultelse	if isstr(thick),		thick= str2num(thick);	endend        [mu,rho,zed,fp,fpp,ray,comp,tran,amu,abs_edge,element_symbol,ifail]=abscal(e_name,energy,thick);if nargout==0,	if zed > 0		s=sprintf(', Z= %2g  Amu= %8.3g  Rho= %8.3g gms/cc',zed,amu,rho);		disp(['Element ',element_symbol,s]);		disp(sprintf ...		('f''= %6.3g f"= %6.3g',fp,fpp));		disp(sprintf ...		('Rayleigh= %6.3g Compton= %8.3g',ray,comp));		disp(sprintf ...		('Absorption length= %6.3g microns at %6.0f eV',1/mu, energy));		disp(sprintf ...		('Transmitted Intensity= %6.3g for %6.3g microns',tran,thick));	else		mu_t= mu* thick;		tran= exp(-mu_t); 		s=sprintf(' has a density of %8.3g gms/cc',rho);		disp([e_name,s]);		disp(sprintf ...		('Absorption length= %6.3g microns at %6.0f eV',1/mu, energy));		disp(sprintf ...		('Transmitted Intensity= %6.3g for %6.3g microns',tran,thick));	endendif nargout> 0,	out_mu= mu;endif nargout> 1,	out_rho= rho;endif nargout== 3,	out_fail= ifail;end
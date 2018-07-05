function [fp,fpp]=henke(zed,energy);% [fp,fpp]=henke(zed[int],energy[eV]);% Get fp, fpp from zed and energy (eV)% based on tables of Henke & Gullickson% Data stored as log10(energy); fp; log10(fpp)% SMB 31-OCT-95global DISK_PATH DISK_DELIM DISK_END ANOMAL_FILES H_LOADED H_ERG H_FP H_FPPenergy=energy(:)';[zed,sym]=element(zed);if min(energy) < 10 | max(energy) > 30000, 	error('Energy out of range'); endif ~isempty( H_LOADED),    idz=find(H_LOADED==zed);else    idz=[];endif isempty( idz), 	idz= length(H_LOADED)+ 1;	H_LOADED(idz)= zed;	load([DISK_PATH,ANOMAL_FILES,DISK_DELIM,'henk_dir',DISK_END,'henke_',sym]);	eval(['H_ERG(idz,:)= ',sym,'_erg;']);	eval(['H_FP(idz,:)= ',sym,'_fp;']);	eval(['H_FPP(idz,:)= ',sym,'_fpp;']);	eval(['clear ',sym,'_erg ',sym,'_fp ',sym,'_fpp']);endlog_erg= log10(energy);for j=1:length(energy)	idx=find(H_ERG(idz,:) <= log_erg(j));	id_lo(j)= idx(length(idx));	idx=find(H_ERG(idz,:) >= log_erg(j));	id_hi(j)= idx(1);endidx=find(id_lo~=id_hi);fp= H_FP(idz,id_lo);fpp= 10.^H_FPP(idz,id_lo);idx_lo=id_lo(idx);idx_hi=id_hi(idx);fract= (log_erg(idx)- H_ERG(idz,idx_lo))./ ...	(H_ERG(idz,idx_hi)- H_ERG(idz,idx_lo));  fp(idx)= H_FP(idz,idx_lo)+ ...	fract.*(H_FP(idz,idx_hi)- H_FP(idz,idx_lo));log_fpp= H_FPP(idz,idx_lo)+ ...	fract.*(H_FPP(idz,idx_hi)- H_FPP(idz,idx_lo));fpp(idx)= 10.^log_fpp;
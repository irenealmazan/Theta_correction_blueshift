%%%  T. Schulli April 2004
%% Petr Mikulik, 20. 4. 2001
%%  .... and Rainer T. Lechner 7.6.2001
%%
%% Runs the unspec program to return h,k,l monitor intensity etc.
%% from spec_file's scan scan_nb. 
%
%%
%% Usage:      s = unspecbm32anghkl( spec_file, scan_nb );
%% Result in:  s.th, s.del, s.gam, s.corrint, s.h, s.k, s.int and TH and INT etc
%%
%%%% the result is saved to a seperate ascii-file
%% 
%% Example:    s15= unspecbm32hkl('suvg.14Apr04, 15 )
%%

function [s, norm] = unspecbm32hkl( spec_file, scan_nb)

if (nargin~=2)
    error('two parameters required');
end

viafajl=1; % read the scan via a temporary file or from stdin
if (viafajl)
  viaf='unspec_tmp -2';
else
  viaf='-';
end
cmd = sprintf('unspec0 %s %s -#  -s  %i -w Eta,Correct,ion2', spec_file, viaf, scan_nb);
%cmd = sprintf('unspec0 %s %s -#  -s  %i -w H,roi1', spec_file, viaf, scan_nb);
%cmd = sprintf('unspec0 %s %s -#  -s  %i -w Energy,detc,attn,Monitor,Detector,Seconds ', spec_file, viaf, scan_nb);
[flag, unspec_tmp] = unix( cmd );
if (flag~=0)
   error( ['COMMAND: ', cmd] );
end
if (viafajl)
  load 'unspec_tmp';
  unix( 'del unspec_tmp' );
end

[nr, nc] = size( unspec_tmp );
s.eta = unspec_tmp(:,1);
s.correct = unspec_tmp(:,2);
s.mon = unspec_tmp(:,3);
%s.detc= unspec_tmp(:,4);
%s.mon=unspec_tmp(:,3);
s.norm=s.correct./s.mon;
%H=s.h; L=s.l; INT=(s.corrint./s.mon2)*(sum(s.mon2)/length(s.mon2));


ausg(:,1)=s.eta;
ausg(:,2)=s.correct;
ausg(:,3)=s.mon;
ausg(:,4)=s.norm;


fid=fopen(['scan' int2str(scan_nb) '.dat'],'w');
%fprintf(fid,'%s \r\n',' H K L I0')
    fprintf(fid,'%2.5e %2.5e 2.5e %2.5e\r\n',ausg');
fclose(fid);




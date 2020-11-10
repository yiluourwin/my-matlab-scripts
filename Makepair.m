% ----------------------------------------------------------------------
% Make pairs for Two-Station Method
% 
% Author: Yi LUO
% Built: 2019-07-12 14:00
% Last Modified
% ----------------------------------------------------------------------

% Input Output
sacpath='sac190711';
pairpath='pairs';

% Parameters
stdist=100; % Minimum distance(km) between two stations 
azdiff=3; % Maximum azimuth for big arc limit

% ----------------------------------------------------------------------

if ~exist(pairpath,'dir')
    mkdir(pairpath)
else
    error('Pairpath exists!')
end

cd(sacpath)
eventlist=dir('./');

for i_event=3:length(eventlist)
    cd(eventlist(i_event).name)
    
    recordlist=dir('./');
    for i_record=3:(length(recordlist)-1)
        sac1=recordlist(i_record).name;
        [sachd1,~]=rsac(sac1);
        for i_record2=(i_record+1):length(recordlist)
            sac2=recordlist(i_record2).name;
            [sachd2,~]=rsac(sac2);
            
            % Distance limit between two stations
            if abs(sachd1.dist-sachd2.dist)<stdist
                continue
            end
            
            % Big arc limit
            A=deg2rad(abs(sachd1.az-sachd2.az));%��������ζ���?
            dist1=deg2rad(sachd1.gcarc);%�������ߵ�Բ�Ľ�
            dist2=deg2rad(sachd2.gcarc);
            
            [B1,B2]=soltri(A,dist1,dist2);%�������?
            
            B1=abs(rad2deg(B1));%�ǻ򲹽�
            B2=abs(rad2deg(B2));
            B=min(B1,B2);
            
            if B>azdiff %�������Բ��ǰ��?
                continue
            end
            
            if exist(['../../',pairpath,'/',strtrim(sachd1.kstnm'),'.',strtrim(sachd2.kstnm')],'dir')
                targetpath=['../../',pairpath,'/',strtrim(sachd1.kstnm'),'.',strtrim(sachd2.kstnm')];
            elseif exist(['../../',pairpath,'/',strtrim(sachd2.kstnm'),'.',strtrim(sachd1.kstnm')],'dir')
                targetpath=['../../',pairpath,'/',strtrim(sachd2.kstnm'),'.',strtrim(sachd1.kstnm')];
            else
                mkdir(['../../',pairpath,'/',strtrim(sachd1.kstnm'),'.',strtrim(sachd2.kstnm')])
                targetpath=['../../',pairpath,'/',strtrim(sachd1.kstnm'),'.',strtrim(sachd2.kstnm')];
            end
            
            copyfile(sac1,targetpath)
            copyfile(sac2,targetpath)
            pairlistID=fopen([targetpath,'/','pairlist'],'a');
            fprintf(pairlistID,'%s\n%s\n',sac1,sac2);
            fclose(pairlistID);
        end
    end
    cd ..
end
cd ..
                
function [B,C]=soltri(A,b,c)
%��������Σ���֪������Բ�Ľ�����нǣ���ʣ�������н�
%��Ϊ����
%
B=atan((sin(A)*sin(b))/(cos(b)*sin(c)-sin(b)*cos(c)*cos(A)));
C=atan((sin(A)*sin(c))/(cos(c)*sin(b)-sin(c)*cos(b)*cos(A)));

end            

close all

%https://www.rohde-schwarz.com/webhelp/ZNB_ZNBT_HTML_UserManual_en/Content/c0e6efab4c3f4280.htm
centFreq = str2double(query(VNA,'FREQ:CENT?'));
fstart = str2double(query(VNA,'FREQ:STAR?'));
fend = str2double(query(VNA,'FREQ:STOP?'));
span = fend - fstart;

%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%

cd 'C:\NEMSMeasure'
tic;

c=fix(clock);
year=c(1); month=c(2); day=c(3); hour=c(4); minute=c(5); seconds=c(6);

%%%%%%%%%%%%%%%%%%%%%%% FILE MANAGEMENT %%%%%%%%%%%%%%%%%%%%%%%
%make data path + copy to shared file
filename = [num2str(year),'-',num2str(month,'%02d'),'-',num2str(day,'%02d'),'-',num2str(hour,'%02d'),num2str(minute,'%02d'),num2str(seconds,'%02d')];
dirname = [num2str(year), '-', num2str(month,'%02d'), '-', num2str(day,'%02d')];
folder = ['C:\matlabdata\' num2str(year) '\' num2str(month,'%02d'),'\JT_auto'];
[status,msg] = mkdir(folder, dirname);
path = [folder, '\', dirname];
copyfile([mfilename, '.m'], [path,'\',filename,'_',mfilename, '.m']); %MFILENAME Name of currently executing MATLAB code file.

%make data path + copy to personnal file
dirname_pers = [num2str(year), '-', num2str(month,'%02d'), '-', num2str(day,'%02d')];
folder_pers = ['\\TW-PHYS.org.aalto.fi\PROJECT\nems\Jonas\matlabdata\JT_auto\' num2str(year) '\' num2str(month,'%02d')];
[status,msg] = mkdir(folder_pers, dirname_pers);
path_pers=[folder_pers, '\', dirname_pers];
copyfile([mfilename, '.m'], [path_pers,'\',filename,'_',extra_name,'_',mfilename, '.m']);

%Making an info file with all the relevant information about the measurement
if exist([path_pers,'\',dirname_pers,'_fileInfo.txt'])==2
  fileInfo = fopen([path_pers,'\',dirname_pers,'_fileInfo.txt'],'a'); % open exist file and append contents
else
  fileInfo = fopen([path_pers,'\',dirname_pers,'_fileInfo.txt'],'w'); % create file and write to it
  fprintf(fileInfo,'File name, Center Freq. (GHz), Start Freq. (GHz), End Freq. (GHz), Span (MHz), Bandwidth (Hz), Power (dBm), Nb of points\n');
end
fprintf(fileInfo,[filename,'_',extra_name,',',num2str(centFreq/1e9),',',num2str(fstart/1e9),',',num2str(fend/1e9),',',num2str(span/1e6),',',num2str(BW),',',num2str(power),',',num2str(numpoints),'\n']);
fclose(fileInfo);
%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%

common.clock=clock;

vari = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

powers = [power];
% powers = linspace(-40,-15,6);

averages(1:length(powers))= 1;

dfreq = linspace(fstart, fend, numpoints);
LinearScale= 0;

TimePerTrace = 2;

timeofmeas = TimePerTrace* length(powers)*mean(averages)/60/60
filename
PlotSpace = 0;  % how many dB offset raw data curves vertically in plots
RunIndex = 1;

MeasTemp = 0;

for pwrind = 1:length(powers)
    if(MeasTemp)
        TNow6 =ResBridgeReadBFchannel(0,6,1);
    else
        TNow6 = 666;
    end
    
    TKelvin6(pwrind) = TNow6;
    
    JT_get_timetrace_ZVB8c;

    ddatamag=magnitude;
    ddata=phase;

    % convert amplitude to lin scale, transmission often looks better
    ddatamaglin = 10.^(ddatamag/20);
    if(LinearScale ~= 0)
        ddatamag = 10.^(ddatamag/20);
    end
    
    % glitches sometimes at first point:
    ddatamag(1)=ddatamag(2);
    ddata(1)=ddata(2);
    
    data(pwrind,:)=ddata;
    datamag(pwrind,:)=ddatamag;
    RunIndex = RunIndex +1;
    
    xaxis = dfreq;
    
    f1 = figure(1);
    subplot(211);
        plot(xaxis, RunIndex*PlotSpace + squeeze(data(pwrind, :)), vari(mod(pwrind,7)+1), 'LineWidth',1);    
    title([filename, ',  uwpwr=', num2str(powers(pwrind))]);
    ylabel('phase'); xlabel('freq'); grid on; hold on; legend(num2str(powers(1:pwrind)'))

    subplot(212);
    plot(xaxis, RunIndex*PlotSpace + squeeze(datamag(pwrind, :)), vari(mod(pwrind,7)+1), 'LineWidth',1);
    title([filename, ',  uwpwr=', num2str(powers(pwrind))]);
    ylabel('mag'); xlabel('freq'); grid on; hold on;legend(num2str(powers(1:pwrind)'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    saveas(f1,[path,'\',filename,'_',extra_name],'fig');
    saveas(f1,[path_pers,'\',filename,'_',extra_name],'fig');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fit S11 dip?
    % Note: Make sure length/delay offset on VNA is set, so phase is flat outside resonance
%     if doFit
%         f2 = figure(100+mod(pwrind,10)); clf;
%         [Qint(pwrind),Qext(pwrind),f0(pwrind)]=S11CircFit2(ddatamaglin.*exp(1i*ddata/180*pi),dfreq,true);
%     end

    %     if powers(pwrind)==10
    pause(1)
    %     end
end

common.clocktotmeastime=toc;
save([path,'\',filename]);%
%%%%%%%%%%%%%%%%%%%%
save([path_pers,'\',filename,'_',extra_name]);
%%%%%%%%%%%%%%%%%%%%
% resetgpib
cd(path);

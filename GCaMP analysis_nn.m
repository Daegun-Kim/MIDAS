clear all;
close all;

Dir = pwd;
ascfile=dir([Dir,'\*.asc']);
numFile=length(ascfile);
sampling_rate=10;
numdatalength=6144;%% 3 blcok = 6144, 9 block = 18432

%% GCaMP ascii file import
for z=1:numFile
    filename=ascfile(z,1).name;
    fid=fopen(filename);
    for i=1:10
        buffer=fgetl(fid);
    end
    raw_data=textscan(fid, '%d %d');
    photon(:,1)=raw_data{1};
    numblock=length(photon)/1024;
    for n=1:numblock
        block(:,n)=photon(1024*(n-1)+1:1024*n,1); % save each block
        F(n,z)=sum(block(:,n));
    end
    F0(:,z)=mean(mean(F(:,z)));   
end
for ii=1:numFile
    total_photon(numblock*(ii-1)+1:numblock*ii,1)=double(F(:,ii));
    result(numblock*(ii-1)+1:numblock*ii,1)=(double(F(:,ii))-F0(:,ii))/F0(:,ii);
end

%% TSPCS test
photon_count_time=linspace(0, 16, 1024);
photon_count_time=photon_count_time';
photon_count(:,1)=mean(block,2);

%% smoothing
GCaMP_smooth=smooth(result, 5);
result = GCaMP_smooth;

%% transient search_GCaMP activity
baseline_window=100;
for iii=1:length(result)-baseline_window
    mu=mean(result(iii:iii+baseline_window-1, 1));
    sigma=std(result(iii:iii+baseline_window-1,1));
    istransient=result(iii+baseline_window,1)-(mu+3*sigma);
    if istransient>=0
        result(iii+baseline_window,2)=true;
    end
end
[r,c]=find(result(:,2)==true);

%% spike sorting of GCaMP activity (result)
spike_times=[r];
diff_spike=spike_times(2:end, 1)-spike_times(1:end-1,1);
sorting_address=find(diff_spike>5); % end of spikes 
start_of_spikes=spike_times(sorting_address+1);
start_of_spikes(2:length(start_of_spikes)+1,1)=start_of_spikes;
start_of_spikes(1,1)=min(spike_times);
end_of_spikes=spike_times(sorting_address);

for p=1:length(end_of_spikes)
    peak(p,1)=max(result(start_of_spikes(p,1):end_of_spikes(p,1),1)); 
    peak_times(p,1)=find(result(start_of_spikes(p,1):end_of_spikes(p,1),1)==peak(p,1))+start_of_spikes(p,1)-1;
end

%% import behavior data
data=xlsread('behavior.xlsx');

sniffing=data(:,1);
sniffing(isnan(sniffing))=[];

%% aligned by sniffing onset
corr_window=40;
for ii=1:length(sniffing(:,1))
    
    single_trace_sniffing(:,ii)=result(sniffing(ii,1)*10-corr_window:sniffing(ii,1)*10+corr_window,1);
    single_time=linspace(-corr_window, corr_window, length(single_trace_sniffing));

end
average_trace_sniffing=mean(single_trace_sniffing,2);
std_trace_sniffing=std(single_trace_sniffing,[],2)/sqrt(length(single_trace_sniffing(1,:)));
figure();
hold on
plot(single_time, average_trace_sniffing, 'r');
hold on
plot(single_time, average_trace_sniffing-std_trace_sniffing, 'k');
hold on
plot(single_time, average_trace_sniffing+std_trace_sniffing, 'k');

sniffing_data(:,1)=average_trace_sniffing;
sniffing_data(:,2)=std_trace_sniffing;

figure();
imagesc(single_trace_sniffing');
color=hot(200);
colormap(color);
caxis([0, 0.05]);

% xlswrite('correlation_trace_data.xlsx', single_trace_sniffing, 'sniffing_raw');
% xlswrite('correlation_trace_data.xlsx', sniffing_data, 'sniffing_avg_std');

%% probability
for i=1:length(single_trace_sniffing(1,:))
    avg_baseline_GCaMP=mean(single_trace_sniffing(1:40,i));
    std_baseline_GCaMP=std(single_trace_sniffing(1:40,i));
    max_GCaMP_sniffing(i,1)=max(single_trace_sniffing(41:end,i));
    if max_GCaMP_sniffing(i,1)>=avg_baseline_GCaMP+2*std_baseline_GCaMP
        success_sniffing(i,1)=1;
    else
        success_sniffing(i,1)=0;
    end
end
probability_sniffing=sum(success_sniffing)/length(success_sniffing);

% % xlswrite('signal_dampening_and_probability.xlsx', max_GCaMP_sniffing,'dampening');
% xlswrite('signal_dampening_and_probability.xlsx', probability_sniffing,'probability', 'B1');
% xlswrite('signal_dampening_and_probability.xlsx', success_sniffing,'probability', 'A1');
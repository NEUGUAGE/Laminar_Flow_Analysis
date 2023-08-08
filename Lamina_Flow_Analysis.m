%%%%%%%%%%%%%%%%%%%%%%Load Data Section Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define ROI
clear data_unique;

%%
%Import the data 
data = readmatrix('choice_20230807_114411_N3_trackedobjects.csv');
%%
%Import the stimulus
stimulus=readmatrix('choice_20230807_114411_stimuli.csv');
%%
%trimming the data using the ROI
x_position_raw=data(:,2);
y_position_raw=data(:,3);
fly_num=3;
if fly_num==1
y_position_raw(y_position_raw<=239-4 | y_position_raw>=827-4)=nan;
x_position_raw(x_position_raw<=310-4 | x_position_raw>=490-4)=nan;

elseif fly_num==2
y_position_raw(y_position_raw<=243-4 | y_position_raw>=819-4)=nan;
x_position_raw(x_position_raw<=523-4 | x_position_raw>=700-4)=nan;

elseif fly_num==3
y_position_raw(y_position_raw<=244-4 | y_position_raw>=820-4)=nan;
x_position_raw(x_position_raw<=734-4 | x_position_raw>=910-4)=nan;

elseif fly_num==4
y_position_raw(y_position_raw<=245-4 | y_position_raw>=822-4)=nan;
x_position_raw(x_position_raw<=944-4 | x_position_raw>=1112-4)=nan;
end
data(:,2)=x_position_raw;
data(:,3)=y_position_raw;


time=data(:,1);
x_position=data(:,2);
y_position=data(:,3);

%%
%find the jumping dot
for i=1:10
if sum(x_position==mode(x_position))>200
    x_position(x_position==mode(x_position))=nan;
end
if sum(y_position==mode(y_position))>200
    y_position(y_position==mode(y_position))=nan;
end
end
for i=1:length(x_position)-1
    if(x_position(i,1)>x_position(i+1,1)*2|x_position(i,1)<x_position(i+1,1)*0.1)
        x_position(i,1)=nan;
    end
    if(y_position(i,1)>y_position(i+1,1)*2|y_position(i,1)<y_position(i+1,1)*0.1)
        y_position(i,1)=nan;
    end
end
%%%%%%%%%%%%%%%%%%%%%%Filter Section Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Replace the NaN Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the indices of the NaN values in the x_position and y_position vectors
nan_indices_x = find(isnan(x_position));
nan_indices_y = find(isnan(y_position));

% Sort the nan_indices_x and nan_indices_y vectors in ascending order
nan_indices_x = sort(nan_indices_x);
nan_indices_y = sort(nan_indices_y);

% Loop over each NaN value and replace it with a linearly interpolated value
for i = 2:length(nan_indices_x)
    % Find the indices of the adjacent non-NaN values and average them
    % replace the NaN
    x_left_index = find(~isnan(x_position(1:nan_indices_x(i))), 1, 'last');
    x_right_index = find(~isnan(x_position(nan_indices_x(i):end)), 1, 'first') + nan_indices_x(i) - 1;
    x_position_index=nan_indices_x(i);
    x_position(x_position_index)=x_position(x_left_index);
end
%Do the same thing for y
for i = 2:length(nan_indices_y)
    % Find the indices of the adjacent non-NaN values
    y_left_index = find(~isnan(y_position(1:nan_indices_y(i))), 1, 'last');
    y_right_index = find(~isnan(y_position(nan_indices_y(i):end)), 1, 'first') + nan_indices_y(i) - 1;
    y_position_index=nan_indices_y(i);
    y_position(y_position_index)=y_position(y_left_index);

end
%%
%replace the data back with x_position
data_unique(:,2)=x_position;
data_unique(:,3)=y_position;
data_unique(:,1)=time;
%%
%replace the first jumping point with the a previous valid data
[~, idx] = unique(data_unique(:,1),'first');
data_unique = data_unique(idx,:);
x_position=data_unique(:,2);
y_position=data_unique(:,3);
time=data_unique(:,1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the velocity
velocity = (sqrt(diff(x_position).*diff(x_position)+diff(y_position).*diff(y_position))./diff(time))/4;%unit conversion 1mm=4px
velocity(isinf(velocity))=nan;
upwind_velocity=(-diff(y_position))./diff(time);%unit conversion 1mm=4px the y axis is inverted
upwind_velocity=upwind_velocity/4;
crosswind_velocity=diff(x_position)./diff(time);
%%
%filter the velocity that is inf and too high indicating jumping dot with NaN
velocity(velocity>20*nanmean(velocity));
upwind_velocity(upwind_velocity>100*nanmean(upwind_velocity));
nan_indices_v = find(isnan(velocity));
nan_indices_v = sort(nan_indices_v);
%%
%Replace the NaN velocity with previous value
for i = 2:length(nan_indices_v)
    % Find the indices of the adjacent non-NaN values
    v_left_index = find(~isnan(velocity(1:nan_indices_v(i))), 1, 'last');
    v_right_index = find(~isnan(velocity(nan_indices_v(i):end)), 1);
    v_position_index=nan_indices_v(i);
    velocity(v_position_index)=(velocity(v_left_index)+velocity(v_right_index))/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%Stimulus%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Matching the stimulus with the time
stimulus_num=size(stimulus,1);%size(frame) of the stimulus
time_num=size(time,1);%size(frame) of the total time
time_stimulus=zeros(size(velocity,1),1);%initialize the time
for i =1:stimulus_num
   for j=1:time_num
        if(time(j,1)>=stimulus(i,1)&&time(j,1)<=stimulus(i,2))
            time_stimulus(j,1)=1;%time that within the stimulus become 1  
        end
   end
end
%Output the time_stimulus with stimulus labeled 1 and non-stimulus with 0
%%
%%%%%%%%%%%%%%%%%%%1.Compute the individual speed(scalar)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% find the average of stimulus velocity and non-stimulus velocity
time_stimulus(1:81,1)=0;
stimulus_velocity=time_stimulus .* velocity;
non_stimulus_velocity=not(time_stimulus).*velocity;
stimulus_upwind_velocity=time_stimulus .* upwind_velocity;
non_stimulus_upwind_velocity=not(time_stimulus).*upwind_velocity;
% calculate the stimulus velocity and non_stimulus velocity
odor_off_time=sum(time_stimulus==0);
odor_on_time=time_num-odor_off_time;%find the stimulus on and off time

odor_v=nanmean(stimulus_velocity);
non_odor_v=nanmean(non_stimulus_velocity);
%%
%column one for speed Δd/Δt
%column two for upwindspeed Δy/Δt
%column three for crosswindspeed Δx/Δt
%column four for stimulus 
output_matrix=zeros(length(velocity),4);
output_matrix(:,1)=velocity;
output_matrix(:,2)=upwind_velocity;
output_matrix(:,3)=crosswind_velocity;
output_matrix(:,4)=time_stimulus;

%%%%%%%%%%%%%%%%%%%%%%Filter Section End%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%Segment Section Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%segment the different trials and select a criteria to
%get rid of the certain inactive trials

%segment the output into different trials
% convert between real time and time frame here
timeframe_stimulus_on=find(time_stimulus==1);
individual_stimulus_duration=10;
pre_stimulus_duration=20;
post_stimulus_duration=40;

time_dot_on_stimuli = find(abs(diff(timeframe_stimulus_on))~=1);

time_dot_on_stimuli = [1;time_dot_on_stimuli];
%time point that stimulus is turning on
real_time_converter=individual_stimulus_duration/(timeframe_stimulus_on(time_dot_on_stimuli(2,1))-timeframe_stimulus_on(time_dot_on_stimuli(1,1)+1)+1);

during_odor_on=timeframe_stimulus_on(time_dot_on_stimuli)-individual_stimulus_duration/real_time_converter;

before_odor=during_odor_on-pre_stimulus_duration/real_time_converter;

during_odor_off=timeframe_stimulus_on(time_dot_on_stimuli);

after_odor=during_odor_off+(post_stimulus_duration)/real_time_converter;
%%
%create a pid delay for data used to measure the flies response when it actually encounter the odor;it takes 4s from top arrive to the bottom
%first get the y_position
pid_delay=(y_position(before_odor)-235)./(810-235)*7;
pid_delay=pid_delay./real_time_converter;

trial_info_matrix=horzcat(before_odor,during_odor_on,during_odor_on+pid_delay+0,during_odor_on+pid_delay+2/real_time_converter,during_odor_off,after_odor);%columns contain before,during and after odor information
%then let's get the upwind_velocity of before after and during in a array
%that composes individual trail
%%
%create a segment upwind 
for i=1:length(trial_info_matrix)
    upwind_velocity_seg_before(:,i)=upwind_velocity(trial_info_matrix(i,1):trial_info_matrix(i,2)-1);
    upwind_velocity_seg_during(:,i)=upwind_velocity(trial_info_matrix(i,3):trial_info_matrix(i,4)-1);
    upwind_velocity_seg_after(:,i)=upwind_velocity(trial_info_matrix(i,5):trial_info_matrix(i,6)-1);
end
%%
%create the same thing for ground speed

for i=1:length(trial_info_matrix)
    velocity_seg_before(:,i)=velocity(trial_info_matrix(i,1):trial_info_matrix(i,2)-1);
    velocity_seg_during(:,i)=velocity(trial_info_matrix(i,2):trial_info_matrix(i,5)-1);
    velocity_seg_after(:,i)=velocity(trial_info_matrix(i,5):trial_info_matrix(i,6)-1);
end
%%
%create the same thing for x,y coordinate
for i=1:length(trial_info_matrix)
    x_position_seg_before(:,i)=x_position(trial_info_matrix(i,1):trial_info_matrix(i,2)-1);
    x_position_seg_during(:,i)=x_position(trial_info_matrix(i,3):trial_info_matrix(i,4)-1);
    x_position_seg_after(:,i)=x_position(trial_info_matrix(i,5):trial_info_matrix(i,6)-1);
end
for i=1:length(trial_info_matrix)
    y_position_seg_before(:,i)=y_position(trial_info_matrix(i,1):trial_info_matrix(i,2)-1);
    y_position_seg_during(:,i)=y_position(trial_info_matrix(i,3):trial_info_matrix(i,4)-1);
    y_position_seg_after(:,i)=y_position(trial_info_matrix(i,5):trial_info_matrix(i,6)-1);
end

%%
%do a t-test on all trials comparision between before stimulus and after
%stimulus

%select for the top 20% within the stimulus to counter the stimulus late
%arrival due to the position of fruit flies in the chamber
%update:the parameter is set to 0.01,which means no effect
threshold_top_20=prctile(upwind_velocity_seg_during,0.01);
upwind_velocity_seg_during_top_20=upwind_velocity_seg_during;
upwind_velocity_seg_during_top_20(:)=0;
for i=1:length(before_odor)
upwind_velocity_seg_during_top_20(:,i)=(upwind_velocity_seg_during(:,i)>threshold_top_20(1,i))
end
upwind_velocity_seg_during_filtered=upwind_velocity_seg_during.*upwind_velocity_seg_during_top_20;
upwind_velocity_seg_during_filtered(upwind_velocity_seg_during_filtered==0)=nan;
for i=1:length(before_odor)
avg_upwind_velocity_before(i,1)=mean(upwind_velocity_seg_before(:,i));
avg_upwind_velocity_during(i,1)=nanmean(upwind_velocity_seg_during_filtered(:,i));
avg_upwind_velocity_after(i,1)=mean(upwind_velocity_seg_after(:,i));
end
%do the same thing for ground speed

for i=1:length(before_odor)
avg_velocity_before(i,1)=mean(velocity_seg_before(:,i));
avg_velocity_during(i,1)=mean(velocity_seg_during(:,i));
avg_velocity_after(i,1)=mean(velocity_seg_after(:,i));
end
%do the same thing 
%%
%create a filter where speed<1mm trials are filtered out
avg_velocity_total=(avg_velocity_during+avg_velocity_before)/2;
avg_upwind_velocity_after(avg_velocity_total<1)=nan;
avg_upwind_velocity_before(avg_velocity_total<1)=nan;
avg_upwind_velocity_during(avg_velocity_total<1)=nan;
%%
%put avg together and do a t-test
total_avg_upwind_velocity=horzcat(avg_upwind_velocity_before,avg_upwind_velocity_during);
[h,p]=ttest(total_avg_upwind_velocity(:,1),total_avg_upwind_velocity(:,2),'Alpha',0.05);
display(h);
if nanmean(velocity)<1.25
    disp("discard,moving speed");
end
if nansum(sqrt(diff(x_position).^2+diff(y_position).^2))/4<50
    disp("discard,moving distance")
end
if sum(isnan(avg_upwind_velocity_during))>0.8*length(avg_upwind_velocity_during)
    disp("discard,not representative")
end
%also based on the graph plotted,if significant tracking error,discard it
%manually,like a lot of straight lines
if h==0
    disp("Not Significant");
else 
    disp("Significant");
end
disp(nanmean(avg_upwind_velocity_before));
disp(nanmean(avg_upwind_velocity_during));
disp(nanmean(avg_upwind_velocity_after));
%%
%Plot the trajectory

time_stimulus_modified=[time_stimulus;0];
stimulus_x=time_stimulus_modified .* data_unique(:,2);
stimulus_y=time_stimulus_modified .* data_unique(:,3);
non_stimulus_x=not(time_stimulus_modified) .* data_unique(:,2);
non_stimulus_y=not(time_stimulus_modified) .* data_unique(:,3);
non_stimulus_x(non_stimulus_x==0)=nan;
non_stimulus_y(non_stimulus_y==0)=nan;
stimulus_x(stimulus_x==0)=nan;
stimulus_y(stimulus_y==0)=nan;
title('Trajectory Plot')
non_stimulus_upwind_velocity_modified=[non_stimulus_upwind_velocity;0];%for plotting purposes,add one 0 to the end, cause by diff()
stimulus_upwind_velocity_modified=[stimulus_upwind_velocity;0];%for plotting purposes,add one 0 to the end, cause by diff()

hold on

plot(stimulus_x,stimulus_y,'r');
plot(non_stimulus_x,non_stimulus_y,'b');
hold off

%%
%{
draw velocity map
hold on
sz=3;
markerType = 'x';
graph2=scatter(stimulus_x,stimulus_y,sz,stimulus_upwind_velocity_modified,markerType); % Plot scatter plot
graph3=scatter(non_stimulus_x,non_stimulus_y,sz,non_stimulus_upwind_velocity_modified);
colormap(copper); % Change colormap
colorbar;
hold off
%}

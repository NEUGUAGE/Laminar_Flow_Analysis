
%%
%%%%%%%%%%Drawing the before after stimulus for effects%%%%%%%%%%%%%%%%
%test plot the graph
for i=1:20
subplot(4,5,i);
hold on
scatter(x_position_seg_before(:,i),y_position_seg_before(:,i),'g','filled');
scatter(x_position_seg_during(:,i),y_position_seg_during(:,i),'r','filled');
scatter(x_position_seg_after(:,i),y_position_seg_after(:,i),'b','filled');
scatter(x_position_seg_during(:,i),y_position_seg_during(:,i),'r','filled');
set(gca,'YDir','reverse')
title(['Trial ' num2str(i)])
hold off
end
%%
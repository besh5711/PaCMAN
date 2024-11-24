function plotESWA(ESWA,eswa,np)
%PLOTESWA Plots eswa and ESWA, or the amplitudes of the exit surface wave
%at the sample plane and detector plane respectively.
%==========================================================================

% Plot ESWA on left side (use to check sampling and edge artifacts):
subplot(1,2,1)
imagesc(ESWA)
title('ESWA')
axis square
set(gca,'fontsize',24)

% Plot eswa on right side (use to check if we have enough scan positions):
subplot(1,2,2)
imagesc(eswa)
title('eswa')
axis square
set(gca,'fontsize',24)
sgtitle(['Position ' num2str(np)],'fontsize',24)
drawnow
end
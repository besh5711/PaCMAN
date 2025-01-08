function plotIterations(fig,NIT,n,obj,prb,P)
%ITERATIONPLOT Makes a plot of the modulus and phase of the object and
%probe, as well as their respective errors for the current iteration.
%==========================================================================
exa_obj = P.exa_obj;
exa_prb = P.exa_prb;

figure(fig);
sgtitle(['Iteration ' num2str(n) ' of ' num2str(NIT)],'fontsize',24)
set(gca,'fontsize',18)

lim1_lb = min(abs(exa_obj),[],'all');
lim1_ub = max(abs(exa_obj),[],'all');

lim2_lb = min(abs(exa_prb),[],'all');
lim2_ub = max(abs(exa_prb),[],'all');

crop = P.cen;

% TOP ROW (Reconstructions) ===============================================
subplot(2,2,1)
imagesc(abs(obj(crop,crop))),title('Rec Obj Modulus','fontsize',18)
colorbar
axis square; axis off
% colormap turbo
% caxis([lim1_lb lim1_ub])

subplot(2,2,2)
imagesc(abs(prb)),title('Rec Prb Modulus','fontsize',18)
% colormap turbo
axis square; axis off
colorbar
% caxis([lim2_lb lim2_ub])

% BOTTOM ROW (Exact) ======================================================
subplot(2,2,3)
imagesc(abs(exa_obj(crop,crop))),title('Exa Obj Modulus','fontsize',18)
colorbar
axis square; axis off
% colormap turbo
caxis([lim1_lb lim1_ub])

subplot(2,2,4)
imagesc(abs(exa_prb)),title('Exa Prb Modulus','fontsize',18)
% colormap turbo
axis square; axis off
colorbar
caxis([lim2_lb lim2_ub])

drawnow;
end
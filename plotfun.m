function p = plotfun(x,y)
figure, 
p = plot(x,y);

% figure, histogram(L8_TOArad); 
xticks([50 100 150 200 250 300 350])
title('Radiance Distribution')
xlabel('Radiance (W/Sr/m^2/{\mum})')
ylabel('Number of Observation')
grid minor
ax  = gca;
ax.FontSize = 36;
end

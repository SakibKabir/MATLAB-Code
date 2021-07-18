rng('default') % For reproducibility
r = normrnd(0, 1, [1 117]);

% looking at the histogram
figure, hist(r)
plot(r, 'o')

RMSE = sqrt(sum(r.^2)/117);
%%
load('L8andAPICS_Libya4')
diff = (Model_Red_band-L8_Red_band)./Model_Red_band;
rms(diff)
figure, plot(diff, 'o')
figure, plot(date, L8_Red_band, 'ro', 'markers', 10)
hold on
plot(date, Model_Red_band, 'b*', 'markers', 10), ylim([0.32 0.36])
%%
diff2 = (Model_Red_band-1.025*L8_Red_band)./Model_Red_band;
rmse2 = rms(diff2);
figure, plot(diff2, 'o')
%%
figure, plot(date, 1.025*L8_Red_band, 'ro', 'markers', 10)
hold on
plot(date, Model_Red_band, 'b*', 'markers', 10), ylim([0.32 0.36])


%% 
L8_Red_band_sim = L8_Red_band;
for i = 1:117
    if L8_Red_band_sim(i) > mean(L8_Red_band)
        L8_Red_band_sim(i) = 1.03*L8_Red_band(i);
    else
        L8_Red_band_sim(i) = 0.97*L8_Red_band(i);
    end
end

%%
figure, plot(date, L8_Red_band_sim, 'ro', 'markers', 10)
hold on
plot(date, Model_Red_band, 'b*', 'markers', 10), ylim([0.32 0.36])

diff2 = (Model_Red_band-L8_Red_band_sim)./Model_Red_band;
rmse2 = rms(diff2)


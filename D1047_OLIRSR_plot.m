%% RSR Plotting Ploting Function: L8 and Dove-R 1047 RSR with Hyperion
function [y] = D1047_OLIRSR_plot()
    %%
    File_L8 = 'Z:\SpecialNeeds\Sakib\Miscellaneous\Ball_BA_RSR.xlsx';
    Dfile_1047 = 'Z:\ImageDrive\PlanetLabs\Processed\1047\Preview\RSRs\1047.csv';
    RSR_1047 = csvread(Dfile_1047, 1, 0);
    % legend 
    lege_1047_L8 = {'1047 Blue','L8 Blue','1047 Green', 'L8 Green','1047 Red','L8 Red', '1047 NIR', 'L8 NIR'};

    %%% Bands and Colors
    bands ={'Blue', 'Green', 'Red', 'NIR'};
    band_colors={'b','g','r','m'};

    y = figure;
    for band =1: 4
    %%% Dove 1047 RSR
        plot(RSR_1047(:,1), RSR_1047(:,band+1), 'color', band_colors{band}, 'LineStyle', '-', 'LineWidth', 2.5)
        hold on

    %%% Landsat 8 RSR 
        L8_RSR = xlsread(File_L8, band+1);
        plot(L8_RSR(:,1), L8_RSR(:,2), 'color', band_colors{band}, 'LineStyle', '--','LineWidth', 1)
        hold on
    end

    xlim([400 1000]);
    ylim([0 1.1]);

    legend(lege_1047_L8,'FontSize',18, 'Location','northeast');
    title('Dove 1047 and L8 RSR with Desert Hyperspectral Profile');
    ylabel('Relative Response')
    xlabel('Wavelength (nm)')

    hold on; grid on; grid minor; ax  = gca; ax.FontSize = 30; ax.GridColor = 'k';
    %   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
end
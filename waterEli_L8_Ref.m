function [L8_TOAref_w_filter] = waterEli_L8_Ref(base_L8, date, percentage) 
    %% Creating OLI Water Elimination Mask from OLI NIR Band
    band = 5;  
    prc = percentage;

    L8_TOAref = TOAref_cal_L8(base_L8, date, band);
    % figure, imagesc(L8_TOAref); colorbar
    % figure, histogram(L8_TOAref)

    % Landsat 8 water filter from NIR band
    L8_TOAref_prc = L8_TOAref;
    L8_TOAref_prc = L8_TOAref_prc(~isnan(L8_TOAref_prc));

    L8_TOAref_prc_sorted = sort(L8_TOAref_prc);
    ind_prc = round(prc*length(L8_TOAref_prc_sorted)); % getting the index at prc
    Data_prc = L8_TOAref_prc_sorted(find(L8_TOAref_prc_sorted) == ind_prc); % getting data at prc

    % Filtering out the data
    L8_TOAref_filtered = L8_TOAref;
    L8_TOAref_filtered(L8_TOAref_filtered <= Data_prc) = NaN; 
    % figure, imagesc(L8_TOAref_filtered); colorbar

    % Creating water elimination mask
    L8_TOAref_w_filter = L8_TOAref_filtered;
    L8_TOAref_w_filter(~isnan(L8_TOAref_w_filter)) = 1;
    % figure, imagesc(L8_TOAref_w_filter); colorbar
end
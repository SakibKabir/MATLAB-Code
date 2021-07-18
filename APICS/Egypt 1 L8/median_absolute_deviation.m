function [ outliersIndex ] = median_absolute_deviation( reflectance,bands,sensitivity )
%{
 The function takes the whole Reflectance array. 
 Define the sensitvity value for  filtering
 bands defines the number of column you want to use.
 This function returns the outliers Row index present in the bands
%}
if(nargin == 1)
   % defining the default band selection from 1 to 7 
   bands = 7;
   % defining the default sensitivity value of 1
   sensitivity = 1;
   warning('Using default value of band = 7 and sensitivity = 1');
elseif (nargin == 2)
    % defining the sensitivity value of 1 if not defined
    sensitivity = 1;
    warning('Using default value of sigma = 1');
else
    % Don't care of this condition
end

% designing the mean absolute deviation mean method 
meanToaSel = reflectance(:,1:bands);

% taking the standard deviation of the Mean TOA Reflectance
stdOfMeanTOA = std(meanToaSel);

% Mean value of the TOA Reflectance of the Each band
meanValue =  mean(meanToaSel);

% compute the absolute differences
absoluteDeviation = abs(meanToaSel - meanValue);

% selecting the median value of the absoluteDeviation
absMedian = median(absoluteDeviation);


% defining the threshold value for filtering the data
thresholdValue = absMedian*sensitivity;

% checking for the outliers
%outliers = absoluteDeviation >  thresholdValue;
outliers = bsxfun(@gt,absoluteDeviation,thresholdValue);

% ROW OF THE EACH OUTLIER PRESENT 
[row,~] = find(outliers == 1);

% removing the duplicate row index
outliersIndex = unique(row);

end


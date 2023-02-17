# qPCR
This script analyses qPCR data from a QuantStudio device.

## Prepare your data files
Export qPCR data as .xls files

## Prepare your settings file
Download the example settings excel file and complete for your samples (one sample per row). You can edit column titles as the program will prompt you to select them. Columns you will need are:

- filename = name of the .xls file.
- well = well number, in 'A01' format
- sample_id = a unique identifier for each sample.
- group1 = grouping variable used to display data in facetted plots (x-axis)
- group2 = grouping variable used to display data in facetted plots (facet columns)
- group3 = grouping variable used to display data in facetted plots (facet rows)
- exclusion_var = column indicating which samples should be excluded the from the analysis.

## Upload data and run analysis
Upload .xls data file/s (you can upload one or multiple) and the excel settings file/s to this app. Select the appropriate columns when prompted.

You will also be asked to input:

- target_source = whether to take determine which reporters relate to which target for each sample from the qPCR data or the settings file.
- calibration_method = select calibrator method for comparative ct analysis. Options are 'Default', 'Callibrator sample/s', 'Lowest GOI ct', or 'Highest GOI ct'. Default uses mean of callibrator sample/s from settings, or if not available then the sample with lowest GOI ct.
- calibrator_var = the column indicating which samples to use as calibrators. Select NA if using the lowest or highest GOI ct as a calibrator method.
- outlier_threshold = this is used in the formula [abs(value - median(value)) > threshold * sd] to identify outlier datapoints.
- remove_outliers = should the wells with outlier ct values be removed from the analysis.
- out.dir = file path to save results.

## Author
Mike Flower, 17/2/2023

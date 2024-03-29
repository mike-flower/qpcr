# qPCR
This script analyses qPCR data from a QuantStudio device. Using this program, you can perform either 'comparative Ct' or 'standard curve' analyses.


## Analysis methods

### Comparative ct (multiple housekeepers)
For comparative ct analyses, this analysis method permits the use of one or more housekeeping genes. There are several ways to work with multiple housekeeping (reference) genes. One way is to select the single best housekeeping gene. Assuming all housekeeping genes work well, and are not affected by experimental conditions, it is also possible to use all to determine relative expression of the gene of interest. This program allows you to do either.

#### The calculations this analysis method uses
Firstly, Ct stands for the cycle threshold (Ct) of your sample. This is given after the qPCR reaction by the qPCR machine. Simply, it is the cycle number where the fluorescence generated by the PCR product is distinguishable from the background noise.

The relative quantity of each target (genes of interest and housekeeping genes) is calculated:

- dCt = calibrator Ct - sample Ct
- Relative quantity = 2^(dCt)

'Calibrator' samples are used to calculate the 'control mean Ct' (CMCT) for each target. Results at the end are therefore relative to the CMCT, and the mean expression of calibrator samples will be 1.

dCt represents how different the Ct of a sample is from the CMCT.

Because PCR amplification (the difference in signal between cycles) is not linear, you need to introduce an exponential term. The '2' in the above equation is the exponential of amplification. 2 assumes a 100% efficient PCR reaction, in which the amount of product doubles every cycle.

Relative gene expression is then:

- Relative gene expression = (relative GOI quantity) / (geomean of relative housekeeping gene quantity)

The above equation is the part that makes use of multiple housekeeping genes. It determines whether the quantity of your gene of interest (relative to the control sample/s) is higher or lower than your housekeepers (each relative to the control sample).


### Comparative ct (2^-ddct)
This is the delta-delta Ct method, also known as the 2<sup>–∆∆Ct</sup> method. It is intended for use with only one housekeeping gene. The formula has been adapted to accomodate multiple housekeeping genes, but the analysis method above [Comparative ct (multiple housekeepers)], would be better suited.

The ∆ (delta) symbol is a mathematical term used to describe the difference between two numbers.

#### The calculations this analysis method uses
First, the technical replicates of each sample are averaged.

Then calculate the ∆Ct for each individual sample. ∆Ct is the difference in Ct values for your gene of interest (GOI) and your housekeeping gene (HK) for a given sample. This is essentially to normalise the gene of interest to a gene which is not affected by your experiment, hence the term housekeeping gene.

- ∆Ct = Ct (GOI) – Ct (HK)

I have adapted this formula to accomodate multiple housekeeping genes, but as above, this is not the best way to handle them.

- ∆Ct = Ct (GOI) – geomean [ Ct (HK) ]

Next, decide which samples to use as 'calibrators'. All results will be expressed relative to these samples. Typically these are control samples, matched to your treated samples - so we would average the ∆Ct values of the control samples to obtain a '∆Ct control mean'. Another way is to pick the sample with the highest or lowest gene of interest Ct (lowest or highest expression level, respectively)

Next, calculate ∆∆Ct for each sample. ∆∆Ct is the difference between the ∆Ct values of the treated samples and the control samples. 

- ∆∆Ct = ∆Ct (treated sample) – ∆Ct (control mean)

Then calculate relative gene expression.

- Relative expression = 2^-(∆∆Ct)


### Standard curve
A serial dilution of standards as comparators to absolutely quantify targets.

#### The calculations this analysis method uses
Firstly, a standard curve of Ct against standard concentration is prepared using the formula.

- ct ~ log (standard concentration)

This model is then used to predict each sample's concentration from it's Ct.

Them, knowing the concentration and the original sample volume, total DNA yield can be calculated.


## Prepare your data for analysis
### Prepare qPCR data files
Export qPCR data as .xls files.

### Prepare your settings file
Download the example settings excel file and complete for your samples (one sample per row). You can edit column titles as the program will prompt you to select them. Columns you will need are:

- filename = name of the .xls file.
- well = well number, in 'A01' format
- sample_id = a unique identifier for each sample.
- group = column containing grouping variable. You can include up to three grouping columns. DO NOT NAME THEM group1, group2 etc. These names are used by the ttest function. Naming convention group_1, group_2 etc is fine.
- reporter1 = the first reporter dye used in each well, e.g. FAM, VIC etc.
- reporter2 ...
- target1 = the target gene that corresponds to reporter1.
- target2 ...
- calibrator_group = groups to which calibrator samples apply. Blank wells form their own group, and if all are left blank then calibrators will be applied to all samples.
- calibrator = column indicating calibrator samples (any character other than NA indicates a calibrator). These are used as calibrators for all samples in the same calibrator group.
- exclude = column indicating which samples to exclude the from the analysis.



## Run the analysis
Download the 'qpcr.R' script and the 'functions' directory, and place both in the same directory on your local machine. Run the script and follow the prompts.

You will be asked to select qPCR data Excel file/s (you can upload one or multiple) and the settings Excel file (again, you can upload one or multiple).

### General inputs you will be asked for
- analysis_method = 'Comparative ct (multiple housekeepers)', 'Comparative ct (2^-ddct)', or 'Standard curve'.
- unique_id = an identifier for each sample. Biological reps will have different identifiers, technical reps will have the same identifier.
- unique_id_label = label for plotting.
- target_source = match reporters (dyes) to targets (genes) using information in the qPCR data file or the settings file.
- empty_samples = samples you expect to be empty, i.e. to fail qPCR (e.g. negative controls). Used to exclude housekeeping genes if amplification fails in a real sample.
- outlier_threshold = this is used in the formula [abs(value - median(value)) > threshold * sd] to identify outlier datapoints.
- remove_outliers = should the wells with outlier ct values be removed from the analysis.
- exclude = column indicating which samples should be excluded the from the analysis.
- pdf attributes = you can set the height and width of the PDF containing plots.
- out.dir = file path to save results.

### Comparative ct inputs
- calibration_method = Choose a method to select calibrator samples. Options are 'Default', 'Calibrator sample/s', 'Lowest GOI ct', 'Highest GOI ct', or 'Select sample/s'. Default uses mean Ct of calibrator sample/s from Settings file, or if not available then the Ct of the sample with lowest gene of interest Ct in each calibrator group.
- calibrator_var = if you selected 'Default' or 'Calibrator samples', this will prompt you to select the column indicating samples to use as calibrators.
- exclude_failed_hk = Select whether to exclude housekeeping targets (genes) if they don't work for all samples. This is recommended; if housekeepers are kept despite failing in some samples, then geomean(HK) for those samples will be skewed towards those HK targets that worked.

### Standard curve inputs
- standard_var = Column containing the standard concentrations.
- conc_unit = default ng/µL
- souce_vol = Column containing the volume of the sample's source. Used to calculate total DNA yield.
- quanfitiler = Indidate whether this is a Quantifiler experiment, in which case you will be prompted to select targets for the Quality index (Small / Large)
- y_log = Choose whether to log scale the dependent variable on plots (y-axis)

### Statistics variables
The program also performs basic t-tests between groups. You can select up to 3 grouping variables (columns in your settings file), which will be presented on the x-axis, facet rows and facet columns respectively.

- hide_ns = Select whether to hide non-significant p-values on the statistics plots.


## Author
Mike Flower, 17/2/2023

# UAT_EFDC
Uncertainty and sensitivity analysis tool for EFDC model
## Modules
1. Perform Latin HyperCube sampling of the parameters to get the input parameter file
2. Read the sampled input parameters, run the model and save the results.
3. Read the output file and calculate the Nash-Sutcliffe coefficient.
## Quick start
### 1. Sampling

![avatar](https://github.com/jlonghku/UAT_EFDC/blob/master/img/Fig1.png)

Figure. 1

1. Input the case path
2. Selcect the range file of sampling parameters, including the parameter name, prior distribution and range of values, see Figure 2 for the specific file format, the first column is the parameter name, the second column is the prior distribution, the third and fourth columns are the minimum and maximum values of the parameters, separated by TAB keys.

![avatar](https://github.com/jlonghku/UAT_EFDC/blob/master/img/Fig2.png)  

Figure. 2

3. Select the input file of EFDC to be modified
4. Input the number of sampling groups
5. Input the random seeds
6. Click the "Run" button to perform the sampling.

The sampling result is generated in file "Parameters.out" and the modified input files are listed in the "inp" folder.

![avatar](https://github.com/jlonghku/UAT_EFDC/blob/master/img/Fig3.png)

Figure. 3

### 2. Run the model

![avatar](https://github.com/jlonghku/UAT_EFDC/blob/master/img/Fig4.png)

Figure. 4

1. Input the model path, it should include all model input files. 
2. Select or edit the "Getefdc.inp". It describes the data that needs to be output. Please see figure 5.

![avatar](https://github.com/jlonghku/UAT_EFDC/blob/master/img/Fig5.png)

Figure. 5

3. Input the total run times
4. Input the threads to be run
5. Click the "Run Model" button

The model will start running and the result will be stored in "RESULT" folder.

### 3. Calculate the Nash-Sutcliffe coefficient

![avatar](https://github.com/jlonghku/UAT_EFDC/blob/master/img/Fig6.png)

Figure. 6

1. Select the monitoring file. It contains the calibration data at specific locations.The first column is the time and other columns are the monitoring results.
2. Select the types of parameters.
3. Input the start and end numder of result to be calculated.
4. Click the "Run Model" button

The tool will start calculating and the result will be stored in "lhv.out" file.

### 4. Settings
Click "Settings" button can setup the path of EFDC.EXE and GetEFDC.exe.

![avatar](https://github.com/jlonghku/UAT_EFDC/blob/master/img/Fig7.png)

Figure. 7
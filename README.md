# Profile-adaptive-sampling
This is the demo code for paper "Multi-profile monitoring based on adaptive sampling". This demo code demonstrates the procedures of parameter estimation, constructing control limit, and obtaining ARL1 for five different methods introduced in the paper.

To replicate the simulation result, run the int_simulation file. We mainly used MASS and Splines packages, also use source function to access the other R file.

Interpret the result: The output of the simulation will be a data frame with proposed, fpca, fpca with \Delta = 0.01,0.1,0.5 in the columns, row 1-7 correspond to r = 8 and mu-shift of 0.4,0.7,1.0,1.4,2.1,2.8,4.2. Similary, row 8-14, correspond to r = 4, row 15-21 correspond to r = 1. SDRL information can be seen in output1 matrix, with every even row represents standard deviation and odd row are ARL_1.

In the demo code, we used only 20 repletions for computing ARL for demonstration purpose. The estimated running time is roughly half an hour depending on the software and hardware. (In the simulation study in our paper, we did 500 times repetitions taking an entire day). If one wants to increase the number of repetition, it's easy to edit (int_simulation file) line 57-67, line 113-120, line 194-202, by changing both i and length of from oc_run04 trough oc_run42.

Details of files and functions:

Method Folder have 4 other files: 1.mfpca_s.R: include two functions

mfpca_est: take data frame and d0 number of dimension to explain in mfpca, output the estimated parameter for mfpca.
m_score: geneating mfpc score for new profile samples

2.dpca_i.R: include two functions

dpca_est:take data frame and d0 number of dimension to explain in fpca, output the estimated parameter for fpca, spitted by list of individual sensors. dpca_score: get fpc score for new profile samples.

3.Get dpcaAi_delta.R: include functions required for ARL_1 for fpca and fpca \Delta = 0.01,0.1,0.5 For fpca: ARLi_initial:Get running length Get_ARLi:Get running length df_hist:Estimate control limit L cusum_df:helper function for redistribution

For foca with delta: ARLi_initial:Get running length Get_ARLi:Get running length df_hist:Estimate control limit L cusum_df:helper function for redistribution

4.Getting ARL_delta.R: include functions required for ARL_1 for proposed method For proposed: ARL_initial: initiate cusum Get_ARL: Get running length mf_hist: Estimate control limit L cusum_stage: helper function for redistribution

Data Folders have all the simulated datasets (put these data along with all function files in the same runing environment when implementing the demo code):

1.Expon_param.Rdata: Model 2 data for both param estimation and control limits training

2.X_param08.Rdata: Model 1 data for both param estimation and control limits training

Link for data used in case study: https://www.kaggle.com/datasets/podsyp/production-quality?select=data_X.csv

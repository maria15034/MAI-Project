# MAI-Project
This is the repository for my final year project at TCD on battery health estimation. 

Here you will find the following,
1. The data generation code which is based on the [P2D model implementation](https://github.com/DEARLIBS) by Dr. Seongbeom Lee.
2. The code I used to pre-process the [experimental data](https://osf.io/qsabn/?view_only=2a03b6c78ef14922a3e244f3d549de78), and also the synthetic data, which I generated using the mentioned P2D model.
3. The code used for modelling it with SVR. 

More details about the data generation process, data pre-processing and data modelling are all included in my thesis. 
For each section, I will include instructions on how to run it and what files are required. Matlab is required for part 1. Parts 2 and 3 are Jupyter Notebooks, so you could use Google Colab for example to run them. The data is downloadable from https://drive.google.com/drive/folders/1fYPlqLbAnzQalBw9erhgFyJrKpaB-VFF?usp=sharing.

Additionally, some code used for analysis, i.e analysis of the metadata (i.e the values for the 8 parameters decided via PSO) are included in the 'Misc.' folder. I also experimented with other models beside SVR, i.e LSTM and TCN, which will be included in this folder.

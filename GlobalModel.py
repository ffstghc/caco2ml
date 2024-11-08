# %% IMPORT PACKAGES + FUNCTIONS

# Data Import
import csv

# Data Transformation 
import numpy as np
import pandas as pd
from statistics import mean 
from scipy.stats import kendalltau

# Cheminformatics
from rdkit.Chem import PandasTools, Descriptors, rdMolDescriptors, AllChem, RDKFingerprint, MACCSkeys
from rdkit.Chem.MolStandardize import rdMolStandardize

# Machine Learning
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFECV
import lightgbm as lgb


# External Caco-2 Test Set (Ro5, TDC)
from tdc.benchmark_group import admet_group
from tdc import Evaluator

# Calcluate all available RDKit descriptors
def getMolDescriptors(mol, missingVal=None):

    res = {}
    for nm,fn in Descriptors._descList:
        # some of the descriptor fucntions can throw errors if they fail, catch those here:
        try:
            val = fn(mol)
        except:
            # print the error message:
            import traceback
            traceback.print_exc()
            # and set the descriptor value to whatever missingVal is
            val = missingVal
        res[nm] = val
        
    return res

# Calculate 1D descriptors of selected set and split into train/validation set    
def select_cmpnd_set(cmpnd_set, flag):
    
    # Choose y (target) data
    if var_model == 'Regression':
        if var_log == 'log10':
            cmpnd_set.replace([np.inf, -np.inf], np.nan, inplace=True) # replace "-inf" by NaN
            cmpnd_set.dropna(subset=['log10 Papp AB Passive'], inplace=True) # Drop NaN rows
            
            y = cmpnd_set['log10 Papp AB Passive']
            
        elif var_log == 'log10cms':
            cmpnd_set.replace([np.inf, -np.inf], np.nan, inplace=True) # replace "-inf" by NaN
            cmpnd_set.dropna(subset=['log10 cm/s Papp,Passive'], inplace=True) # Drop NaN rows
            
            y = cmpnd_set['log10 cm/s Papp,Passive']
            
        else: # nm/s
            y = cmpnd_set['Papp Passive [nm/s]']
            
    elif var_model == 'Classification':
        
        cmpnd_set['Papp Passive Category'] = cmpnd_set.apply(label_permeability, axis=1) # Label category
        
        if var_log == 'log10':
            y = cmpnd_set['Papp Passive Category']
        else:
            y = cmpnd_set['Papp Passive Category']
        
    # ECFP6 fingerprint (control)
    if flag == 'FP': 
        X, FP = calc_1D_descriptors(cmpnd_set) 
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state=94)
    
        return X_train, X_test, y_train, y_test, X, FP
    
    # RDKit fingerprint
    elif flag == 'FP_RDKit':
        X, FP = calc_1D_descriptors(cmpnd_set) 
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state=94)
    
        return X_train, X_test, y_train, y_test, X, FP   
    
    elif flag == 'MACCS':
        X, FP = calc_1D_descriptors(cmpnd_set) 
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state=94)
    
        return X_train, X_test, y_train, y_test, X, FP  

    # Imported, pre-calcualted PaDEL descriptors
    elif flag == 'PaDEL':
        X = pd.read_csv('') # Local path to precomputed PaDEL descriptors removed

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state=94)
    
        return X_train, X_test, y_train, y_test
    
    # All available RDKit descriptors
    elif flag == '0D_RDKit':
        
        # Calculate RDKit descriptors for full data to allow export and extraction 
        All_0D_RDKit = [getMolDescriptors(m) for m in df_imp_FULL['Molecule STD']]        
        X = pd.DataFrame(All_0D_RDKit) # Full dataframe 
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state=94)
    
        return X_train, X_test, y_train, y_test
    
    elif flag == 'CDDD':
        
        unique_smiles = cmpnd_set['Structure STD']

        cddd_server = CDDDRequest(port=) # Port to Bayer internal CDDD API removed

        X = pd.DataFrame(cddd_server.smiles_to_cddd(list(unique_smiles),preprocess=False)) # Convert smiles to CDDD descriptor

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=94)
        
        return X_train, X_test, y_train, y_test

# Label Papp Passive Category
def label_permeability(row):
    
    if row['Papp Passive [nm/s]'] < 10:
        return 'Low'
    
    elif 10 <= row['Papp Passive [nm/s]'] <= 70:
        return 'Medium'
    
    else:
        return 'High'

# Calculation of 1D fingerprint descriptors
def calc_1D_descriptors(dataframe):
    
    if var_descriptor_set == 'FP':    # ECFP (Circular)
        
        # Parameters
        radius = 1          # Std: 2
        nBits = 1024        # Std: 1024

        # Loop all Molecules
        ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x, radius, nBits) for x in dataframe['Molecule STD']]

        # Create DF with Index and Fingerprints
        ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
        ecfp6_bits = [list(l) for l in ECFP6]
        df_ECFP6 = pd.DataFrame(ecfp6_bits, index = dataframe.index, columns=ecfp6_name)
        df_ECFP6.head(1)   # Show First Entry of Bit Output
            
        return df_ECFP6, ECFP6
    
    elif var_descriptor_set == 'FP_RDKit': # RDKit (Topological)
        
        # Parameters
        length = 5
        nBits = 1024
        
        fpgen = AllChem.GetRDKitFPGenerator(maxPath=length, fpSize=nBits)
        RDKitFP = [fpgen.GetFingerprint(x) for x in dataframe['Molecule STD']]

        RDKitFP_name = [f'Bit_{i}' for i in range(nBits)]
        RDKitFP_bits = [list(l) for l in RDKitFP]
        df_RDKitFP = pd.DataFrame(RDKitFP_bits, index = dataframe.index, columns=RDKitFP_name)
        
        return df_RDKitFP, RDKitFP_bits
    
    
    elif var_descriptor_set == 'MACCS': # MACCS (Structural):
        
        MACCS_Keys = [rdMolDescriptors.GetMACCSKeysFingerprint(x) for x in dataframe['Molecule STD']]
        MACCS_bits = [list(l) for l in MACCS_Keys]
        df_MACCS_Keys = pd.DataFrame(MACCS_bits, index = dataframe.index)
        
        return df_MACCS_Keys, MACCS_bits

# Calculation of 1D fingerprint descriptors for TDC benchmark
def calc_1D_descriptors_Benchmark(dataframe):
    
    # Add moelcules object
    PandasTools.AddMoleculeColumnToFrame(dataframe,'Drug','Molecule',includeFingerprints=True)
    PandasTools.RemoveSaltsFromFrame(dataframe, molCol='Molecule')
    
    if var_Benchmark_Descriptor == "ECFP":
            
        # Parameters
        radius = 3
        nBits = 2048

        # Loop all Molecules
        ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x, radius, nBits) for x in dataframe['Molecule']]

        # Create DF with Index and Fingerprints
        ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
        ecfp6_bits = [list(l) for l in ECFP6]
        df = pd.DataFrame(ecfp6_bits, index = dataframe.index, columns=ecfp6_name)
        
    elif var_Benchmark_Descriptor == "RDKit209":
        
        All_0D_RDKit = [getMolDescriptors(m) for m in dataframe['Molecule']]        
        df = pd.DataFrame(All_0D_RDKit) # Full dataframe 
        
    elif var_Benchmark_Descriptor == "CDDD":
        
        unique_smiles = dataframe['Drug']
        cddd_server = CDDDRequest(port=) # Port for Bayer internal CDDD API removed
        df = pd.DataFrame(cddd_server.smiles_to_cddd(list(unique_smiles),preprocess=False))
        
    return df


# %% CALCULATE SELECTED DESCRIPTORS AND AUTOMATIC TRAIN TEST SPLIT ##
#####################################################################

## Choose descriptor ('FP', 'PhysChem')  
## Choose data transformation ('log10', 'log10cms, '')
## Choose classification/regression y data 
var_df_compounds = df_imp_FULL      	# df_imp_Full is dataframe with Caco-2 data and molecular structures
var_descriptor_set = 'FP'               # FP, PaDEL, FP_RDKit, MACCS, 0D_RDKit, CDDD
var_log = 'log10cms'                    # log10, '', log10cms
var_model = 'Regression'                # Regression, Classification
##########################################################################

if var_descriptor_set == 'FP':
    X_train, X_test, y_train, y_test, df_ECFP6, lst_ECFP6 = select_cmpnd_set(var_df_compounds, 'FP')  
elif var_descriptor_set == 'FP_RDKit':
    X_train, X_test, y_train, y_test, df_RDKIT_FP, lst_RDKIT_FP = select_cmpnd_set(var_df_compounds, 'FP_RDKit')
elif var_descriptor_set == 'MACCS':
    X_train, X_test, y_train, y_test, df_MACCS_keys, lst_MACCS_keys = select_cmpnd_set(var_df_compounds, 'MACCS') 
elif var_descriptor_set == 'CDDD':
    X_train, X_test, y_train, y_test = select_cmpnd_set(var_df_compounds, 'CDDD')
elif var_descriptor_set == 'PaDEL':
    X_train, X_test, y_train, y_test = select_cmpnd_set(var_df_compounds, 'PaDEL')
elif var_descriptor_set == '0D_RDKit':
    X_train, X_test, y_train, y_test = select_cmpnd_set(var_df_compounds, '0D_RDKit')


# %% SVM REGRESSOR ##
#####################

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)


# Configure SVM
# regr_SVM_rbf = SVR(kernel="poly", C=100, gamma="auto", degree=3, epsilon=0.1, coef0=1) # 4 h execution duration

# regr_SVM_rbf = SVR(kernel="rbf") # 0.5 h duration

regr_SVM_rbf = SVR(kernel='rbf') # 

# Fit
regr_SVM_rbf.fit(X_train_scaled, y_train)

# Predict
y_pred = regr_SVM_rbf.predict(X_test_scaled)   


# %% MLP REGRESSOR ##
#####################

# Standard scale data
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Configure MLP
regr_MLP = MLPRegressor(random_state=1, 
                                      max_iter=1000, 
                                      verbose=2, 
                                      hidden_layer_sizes = (1024), # 1024, 600
                                      early_stopping=True)

# Train MLP
regr_MLP.fit(X_train, y_train)

# Predict MLP
y_pred = regr_MLP.predict(X_test)


# %% LIGHTGBM REGRESSION MODEL (INTERNAL TEST SET) ##
####################################################

# Configure LightGBM
params = {
    'boosting': 'gbdt',
    'objective': 'regression',
    'num_leaves': 35,
    'n_estimators' : 2000,
    'learning_rate': 0.05,
    'metric': {'l1','l2'},
    'verbose': -1,
    'n_jobs' : 8
}

# Train (and cross-validate if wanted) LightGBM
lgb_train = lgb.Dataset(X_train, y_train)
# lgb_eval = lgb.Dataset(X_test, y_test, reference=lgb_train)

rgr_LightGBM = lgb.train(params,
                 train_set=lgb_train) # ,valid_sets=lgb_eval)

cv_mod = lgb.cv(params, 
                lgb_train, 
                500, 
                nfold = 10,
                stratified = False)

# Convert to CV result table
df_CV_LightGBM = pd.DataFrame.from_dict(cv_mod, orient='index').T
df_CV_LightGBM = df_CV_LightGBM.rename(columns={"valid l1-mean": "Mean MAE",
                                                "valid l1-stdv": "MAE SD",
                                                "valid l2-mean": "Mean RMSE",
                                                "valid l2-stdv": "RMSE SD"})

# Print CV results
print("Overall Mean MAE after 10-fold CV: ", np.mean(df_CV_LightGBM['Mean MAE']))
print("Overall Mean MAE SD after 10-fold CV: ", np.mean(df_CV_LightGBM['MAE SD']))
print("Overall Mean RMSE after 10-fold CV: ", np.sqrt(np.mean(df_CV_LightGBM['Mean RMSE'])))
print("Overall Mean RMSE SD after 10-fold CV: ", np.sqrt(np.mean(df_CV_LightGBM['RMSE SD'])))

# Prediction for test set
y_pred = rgr_LightGBM.predict(X_test)


# %% TDC EXTERNAL TEST BENCHMARK ##
###################################

##########################
### CHANGE VALUES HERE ###
##########################
var_Benchmark_Regressor     = 'LightGBM'    # MLP LightGBM SVM
var_Benchmark_Train         = 'external'    # external internal
var_Benchmark_Descriptor    = 'RDKit209'    # ECFP RDKit209 CDDD
var_log                     = 'log10cms'    # nm/s log 
###################################################### 
    
# Import TDC data
group = admet_group(path = 'data/')
predictions_list = []

for seed in [1, 2, 3, 4, 5]:
    benchmark = group.get('Caco2_Wang') 
    
    ## all benchmark names in a benchmark group are stored in group.dataset_names
    predictions = {}
    name = benchmark['name']
    train_val, test = benchmark['train_val'], benchmark['test']
    train, valid = group.get_train_valid_split(benchmark = name, split_type = 'default', seed = seed)
    
    Benchmark_TrainVal = calc_1D_descriptors_Benchmark(train_val)
    Benchmark_Test = calc_1D_descriptors_Benchmark(test)
    
    ## Data conversion (log, non-log)
    if var_log == 'nm/s':
        test['Y'] = np.power(10, test['Y'])     # De-log
        test['Y'] = test['Y'].mul(10**7)    # Convert cm/s to nm/s
        train_val['Y'] = np.power(10, train_val['Y'])   # De-log
        train_val['Y'] = train_val['Y'].mul(10**7)  # Convert cm/s to nm/s
        
    elif var_log == 'log': # Only affects internal data
        y_test = y_test * 10**-7   # Convert non-log nm/s to log cm/s
        y_train = y_train * 10**-7  # Convert non-log nm/s to log cm/s
        
    
    
    ## Model selection
    if var_Benchmark_Regressor == 'LightGBM':
    
        # LightGBM Parameters
        params = {
            'boosting': 'gbdt',
            'objective': 'regression',
            'num_leaves': 35,
            'n_estimators' : 2000,
            'learning_rate': 0.05,
            'metric': {'l1','l2'},
            'verbose': -1}

        ## Training data selection
        if var_Benchmark_Train == 'external':
            lgb_train = lgb.Dataset(Benchmark_TrainVal, train_val['Y'])
            lgb_eval = lgb.Dataset(Benchmark_Test, test['Y'], reference=lgb_train)
            
        elif var_Benchmark_Train == 'internal':
            lgb_train = lgb.Dataset(X_train, y_train)
            lgb_eval = lgb.Dataset(X_test, y_test, reference=lgb_train)

        var_regr = lgb.train(params,
                    train_set=lgb_train) # ,valid_sets=lgb_eval)
    
    elif var_Benchmark_Regressor == 'MLP':
        
        from sklearn.neural_network import MLPRegressor

        var_regr = MLPRegressor(random_state=1, max_iter=1000, verbose=2,
                            hidden_layer_sizes = (600, 2048), early_stopping=True)
        
        if var_Benchmark_Train == 'external':
            var_regr.fit(Benchmark_TrainVal, train_val['Y'])
               
        elif var_Benchmark_Train == 'internal':
            var_regr.fit(X_train, y_train)
            
            
    elif var_Benchmark_Regressor == 'SVM':
        
        from sklearn.svm import SVR
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        var_regr = make_pipeline(StandardScaler(),
                     SVR(kernel="linear", C=1, gamma=0.1, epsilon=0.1, verbose=1))
        
        if var_Benchmark_Train == 'external':
            var_regr.fit(Benchmark_TrainVal, train_val['Y'])
               
        elif var_Benchmark_Train == 'internal':
            var_regr.fit(X_train, y_train) 


    ## Predict for Test Set
    y_pred_test = var_regr.predict(Benchmark_Test)
        
    predictions[name] = y_pred_test
    predictions_list.append(predictions)

results = group.evaluate_many(predictions_list)
print(results)

# Evaluation
lst_metrices = ['MSE','MAE', 'RMSE', "R2", "PCC", "Spearman"]


for i in lst_metrices:
    
    evaluator = Evaluator(name = i)
    score = evaluator(test['Y'], y_pred_test)
    
    print(i, score)



# %% RECURSIVE FEATURE ELIMINATION WITH LIGHTGBM ##
###################################################

X_full = pd.concat([X_train, X_test])
y_full = pd.concat([y_train, y_test])

min_features_to_select = 1  # Minimum number of features to consider

# Configure LightGBM for RCFECV
rgr = lgb.LGBMRegressor(boosting_type='gbdt', 
                        objective='regression',
                        num_leaves=35, 
                        n_estimators=2000, 
                        learning_rate=0.05,
                        n_jobs=8) 

# Configure cross-validation for RFE
cv = KFold(5)

# Configure recurssive feature elimination (RFE)
rfecv = RFECV(
    estimator=rgr,
    step=1,
    cv=cv,
    scoring="neg_mean_absolute_error",
    min_features_to_select=min_features_to_select,
    n_jobs=8,
)

# Fit RFE
rfecv.fit(X_full, y_full)

# Output RFE fit results
print(f"Optimal number of features: {rfecv.n_features_}")


# %%  LightGBM Classification Model ##
######################################

# Configure LGBM classifier
clf_LightGBM = lgb.LGBMClassifier(boosting_type= 'gbdt',
                                  device = 'gpu',
                                  objective = 'multiclass',
                                  max_depth= 0,
                                  num_leaves = 35,
                                  learning_rate = 0.05,
                                  n_estimators= 2000,
                                  class_weight= 'balanced',
                                  n_jobs= 8) 

# Fit LGBM classifier
clf_LightGBM.fit(X_train, y_train) 

# Predict with LGBM classifier
y_pred = clf_LightGBM.predict(X_test)


# %% 10-Fold Cross-Validation ##
################################

# Standardscale
'''
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)
'''

# Combine Train + Test for Full Data
X = pd.concat([X_test, X_train])
y = pd.concat([y_test, y_train])

# Reset Index
X = X.reset_index(drop=True)
Y = y.reset_index(drop=True)

# Setup 10-Fold CV mit Mean MAE + SD Calculation
kf = KFold(n_splits=10)

# MLP Parameter
regr_MLP = MLPRegressor(random_state=1,
                        max_iter=1000,
                        verbose=2,
                        hidden_layer_sizes = (600), # Deep: 1024, 600 Shallow: 600
                        early_stopping=True)

# LightGBM Parameter
params = {
    'boosting': 'gbdt',
    'objective': 'regression',
    'num_leaves': 35,
    'n_estimators' : 2000,
    'learning_rate': 0.05,
    'metric': {'l1','l2'},
    'verbose': -1,
    'n_jobs' : 8}

# SVM Parameter
# regr_SVM_rbf = SVR(kernel='rbf')
regr_SVM_poly = SVR(kernel='poly')

lst_MAEs_10CV = []
lst_R2s_10CV = []
lst_Taus_10CV = []


for i, (train_index, test_index) in enumerate(kf.split(X)):
    print(f"Fold {i}:")
    print(f"  Train: index={train_index}")
    print(f"  Test:  index={test_index}")
    
    
    #########################################################
    # REMOVE/ADD COMMENTS BELOW TO SELECT APPLICABLE MODELS #
    #########################################################
    '''
    # Fit MLP ########################################################
    # regr_MLP.fit(X.iloc[train_index], y.iloc[train_index])

    # y_pred = regr_MLP.predict(X.iloc[test_index])
    ##################################################################
    
    
    # Fit LightGBM ###################################################
    lgb_train = lgb.Dataset(X.iloc[train_index], y.iloc[train_index])

    rgr_LightGBM = lgb.train(params,
                    train_set=lgb_train) # ,valid_sets=lgb_eval)

    # Prediction for test set
    y_pred = rgr_LightGBM.predict(X.iloc[test_index])
    ##################################################################
    

    # Fit RBF SVM ####################################################
    regr_SVM_rbf.fit(X.iloc[train_index], y.iloc[train_index])

    # Predict
    y_pred = regr_SVM_rbf.predict(X.iloc[test_index])
    
    '''
    
    # Fit Poly SVM ###################################################
    regr_SVM_poly.fit(X.iloc[train_index], y.iloc[train_index])

    # Predict
    y_pred = regr_SVM_poly.predict(X.iloc[test_index])
    ##################################################################

    # Collect Metrics
    lst_MAEs_10CV.append(round(mean_absolute_error(y.iloc[test_index], y_pred), 3))
    lst_R2s_10CV.append(round(r2_score(y.iloc[test_index], y_pred), 3))
    lst_Taus_10CV.append(kendalltau(y.iloc[test_index], y_pred))


# Extract Kendall-Tau Results
for i in range(0,len(lst_Taus_10CV)):
    lst_Taus_10CV[i] = round(lst_Taus_10CV[i].statistic, 3)    

# Combine collected Metrics    
df_10CV_Metrics = pd.DataFrame(data = {'MAE (10 CV)': lst_MAEs_10CV,
                                       'R2 (10 CV)': lst_R2s_10CV,
                                       'K-Tau (10 CV)': lst_Taus_10CV})
# Calculate Mean and SD for Metrics
df_10CV_Metrics_Means_SDs = {'MAE (Mean 10 CV)': round(mean(lst_MAEs_10CV), 3),
                             'MAE (SD 10 CV)' : round(np.std(lst_MAEs_10CV), 3),
                             'K-Tau (Mean 10 CV)': round(mean(lst_Taus_10CV), 3),
                             'K-Tau (SD 10 CV)': round(np.std(lst_Taus_10CV), 3),
                             'R2 (Mean 10 CV)': round(mean(lst_R2s_10CV), 3),
                             'R2 (SD 10 CV)' : round(np.std(lst_R2s_10CV), 3)}


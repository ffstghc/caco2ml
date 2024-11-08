# %% LOCAL ADAPTIVE MODEL (FIXED SIMILAR STRUCTURES) ##
#######################################################

import math
import pandas as pd
import numpy as np
import lightgbm as lgb 
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# FUNCTIONS

def calc_similarity_pair(a, b):
    
    if a is None or b is None: 
        return 0.0
    amol = Chem.MolFromSmiles(a)
    bmol = Chem.MolFromSmiles(b)
    if amol is None or bmol is None:
        return 0.0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(amol, 3, nBits=2048, useChirality=False)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(bmol, 3, nBits=2048, useChirality=False)
        
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    return similarity

def calc_regression_metrics(y_test_data, y_predicted):
    
    y_test_data = pd.DataFrame(y_test_data)
    y_predicted = pd.DataFrame(y_predicted)

    y_test_data = y_test_data.dropna()
    y_predicted = y_predicted.dropna()
    
    MSE =  round(mean_squared_error(y_test_data, y_predicted), 3)
    MAE = round(mean_absolute_error(y_test_data, y_predicted), 3)
    RMSE = round(math.sqrt(mean_squared_error(y_test_data, y_predicted)), 3)
    R2 = round(r2_score(y_test_data, y_predicted), 3)

    return MSE, MAE, RMSE, R2

def LGBM_Regression(X_train, y_train, X_test):
    # MODEL PARAMETERS
    lgb_train = lgb.Dataset(X_train, y_train)

    rgr_LightGBM = lgb.train(params,
                    train_set=lgb_train)

    # Prediction for test set
    y_pred = rgr_LightGBM.predict(X_test)
    
    return y_pred

# LGBM model paramter
params = {
    'boosting': 'gbdt',
    'objective': 'regression',
    'num_leaves': 35,
    'n_estimators' : 2000,
    'learning_rate': 0.05,
    'metric': {'l1','l2'},
    'verbose': -1,
    'n_jobs' : 16,
    'random_state': 123
    }

# Select fixed cut-offs for similarity to select data to train with
lst_CutOffs = range(0, 0.995, 0.005)

# Loop for all cutoffs
for cutoff in lst_CutOffs:

    # Initialize lists
    lst_Y_Predictions = []    
    lst_Y_Original = []
    lst_amount_training_structures = []
    lst_max_Similarity = []
    lst_mean_Similarity = []

    # Loop for all SMILES

    for i in range(0, len(df_Test)):    #

        # Remove test compound from training
        print('Remove SMILES from training:', df_Test['Structure STD'].iloc[i])  
        print('n(Training) BEFORE test cmpnd removal from training data: ', len(df_Train_Pool))
        
        df_Train_Pool_NO_TEST = df_Train_Pool.drop(df_Train_Pool[df_Train_Pool['Structure STD'].isin([df_Test['Structure STD'].iloc[i]])].index, axis='index')

        print('n(Training) AFTER test cmpnd removal from training data: ', len(df_Train_Pool_NO_TEST)) 

        
        print('Calculating similarities between test compound and pool...')
        df_Train_Pool_NO_TEST['Similarity'] = None

        for row in range(0, len(df_Train_Pool_NO_TEST)):
            df_Train_Pool_NO_TEST['Similarity'].iat[row] = calc_similarity_pair(df_Test['Structure STD'].iloc[i],
                                                                                        df_Train_Pool_NO_TEST['Structure STD'].iloc[row])


        #_________________________________________________________
        ## Find compound cluster with selected tanimoto similarity
        print('Extracting compounds meeting the similarity cut-off:', cutoff, '%')

        # Initialize lists
        lst_sim_ext_inds = []
        lst_sim_ext_sim_values = []

        try:
            for ind in df_Train_Pool_NO_TEST.index.values:
                if df_Train_Pool_NO_TEST['Similarity'].loc[ind] > cutoff:
                    lst_sim_ext_sim_values.append(df_Train_Pool_NO_TEST['Similarity'].loc[ind])    # Actual similarity values
                    lst_sim_ext_inds.append(ind)                            # Index of compounds meeting cutoff (0 - 35k)

            print('Highest similarity: ', max(lst_sim_ext_sim_values))  
        except:
            print('No highest similarity found')  
            lst_sim_ext_sim_values.append(None)   
            lst_sim_ext_inds.append(None)   
            print('Similarities and train indices replaced by "None"')        
    
        
        #_________________________________________________________
        ## DESCRIPTOR CALCULATION + TRAIN

        try:
            # Select training and data
            print("Extracting 0D RDKit training descriptors from pool for test compound ", i)
            X_train = df_Train_Pool_NO_TEST.loc[:, 'MaxAbsEStateIndex':'fr_urea']
            X_train = X_train.loc[lst_sim_ext_inds]   
            y_train = df_Train_Pool_NO_TEST['log10(Papp AB) [cm/s]'].loc[lst_sim_ext_inds] 
        except:
            print('Training descriptor extraction failed for test compound ', i)
            X_train = None
            y_train = None
            print('Train Data replaced by "None"')

        try:    
            print("Extracting 0D RDKit descriptors for test compound ", i)
            X_test = df_Test.loc[:, 'MaxAbsEStateIndex':'fr_urea']
            X_test = X_test.iloc[i]
            y_test = df_Test['log10(Papp AB) [cm/s]'].iloc[i]

        except:
            print('Test descriptor extraction failed for test compound ', i)
            X_test = None
            y_test = df_Test['log10(Papp AB) [cm/s]'].iloc[i]
            print('Test Data replaced by "None"')
            


        #_________________________________________________________
        ## Train model on similar compounds + Predict
        ## LightGBM 
        try:
            print('Predicting for test compound #', i)

            y_pred = LGBM_Regression(X_train, y_train, X_test)
            
            print('Predicted value:', y_pred)
            print('Experimental value:', y_test)
        except:
            print('Predicting failed for test compound', i )
            y_pred = None
            print('y_pred replaced by "None"')


        try:    
            # Collect data
            lst_Y_Predictions.append(y_pred[0])  # y_pred (predicted)
            lst_Y_Original.append(y_test)
        except:    
            print('Appending Pred data failed, filling NA')
            lst_Y_Predictions.append(None)   # y_pred (predicted)  
            lst_Y_Original.append(None)  
                                                                             
               
        try:    
            lst_amount_training_structures.append(len(lst_sim_ext_inds))                # n(Training)
        except:
            print('Appending n(Training) failed, filling NAs')
            lst_amount_training_structures.append(None) 
        

        try:    
            lst_sim_ext_sim_values_notNone = [x for x in lst_sim_ext_sim_values if x is not None]
            lst_max_Similarity.append(max(lst_sim_ext_sim_values_notNone))                      # Maximum similarity in training data
        except:    
            print('Appending similarity data failed, filling NAs')            
            lst_max_Similarity.append(None)                                             # Maximum similarity in training data


        try:    
            # Remove None elements
            lst_sim_ext_sim_values_notNone = [x for x in lst_sim_ext_sim_values if x is not None]
            lst_mean_Similarity.append(np.mean(lst_sim_ext_sim_values_notNone))                 # Mean similarity in training data
        except:
            lst_mean_Similarity.append(None)                                            # Mean similarity in training data
       

        print('\n') # New line    

    try:
    # Combine results per cutoff
        df_Y_test = pd.DataFrame({'Predicted': lst_Y_Predictions, 
                                    'Measured': lst_Y_Original,
                                    'n': lst_amount_training_structures,
                                    'max. Similarity': lst_max_Similarity,
                                    'mean Similarity': lst_mean_Similarity})
    except:
        print('Combining results failed.')
        

    try:    
        # Collect metrics
        var_MSE, var_MAE, var_RMSE, var_R2 = calc_regression_metrics(df_Y_test['Measured'], df_Y_test['Predicted'])
        dict_metrics = {'MSE': var_MSE, 'MAE': var_MAE, 'RMSE': var_RMSE, 'R2': var_R2}
        df_Metrics = pd.DataFrame([dict_metrics])

        print(df_Metrics)
        print(df_Y_test)
    except:
        print('Collecting metrics failed, filling with NAs')
        dict_metrics = {'MSE': None, 'MAE': None, 'RMSE': None, 'R2': None}
        df_Metrics = pd.DataFrame([dict_metrics])


    try:
        # Export metrics
        file_name_Metrics = 'Metrics_' + str(cutoff) + '_0DRDKit_log10cms_LGBM_FixedSimEval'
        df_Metrics.to_pickle(file_name_Metrics)
    except:
        print('Exporting metrics failed.')


    try:
        # Export results
        file_name_Results = 'Results_' + str(cutoff) + '_0DRDKit_log10cms_LGBM_FixedSimEval' # File name dependent on cutoff
        df_Y_test.to_pickle(file_name_Results) # Save predictions in file
    except: 
        print('Exporting results failed.')

    

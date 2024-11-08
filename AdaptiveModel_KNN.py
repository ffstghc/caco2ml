# %% LOCAL ADAPTIVE MODEL (KNN TRAINING DATA SELECTION) ##
##########################################################

import math
import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.neighbors import NearestNeighbors
from scipy.stats import kendalltau



def calc_regression_metrics(y_test_data, y_predicted):
    
    y_test_data = pd.DataFrame(y_test_data)
    y_predicted = pd.DataFrame(y_predicted)
    y_test_data = y_test_data.dropna()
    y_predicted = y_predicted.dropna()
    
    MSE =  round(mean_squared_error(y_test_data, y_predicted), 3)
    MAE = round(mean_absolute_error(y_test_data, y_predicted), 3)
    RMSE = round(math.sqrt(mean_squared_error(y_test_data, y_predicted)), 3)
    R2 = round(r2_score(y_test_data, y_predicted), 3)
    k_tau = kendalltau(y_test_data, y_predicted)

    return MSE, MAE, RMSE, R2, k_tau

def LGBM_Regression(X_train, y_train, X_test):
    # MODEL PARAMETERS
    lgb_train = lgb.Dataset(X_train, y_train)

    rgr_LightGBM = lgb.train(params,
                    train_set=lgb_train)

    # Prediction for test set
    y_pred = rgr_LightGBM.predict(X_test)
    
    return y_pred

# MODEL PARAMETER
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

# Select cut-offs for similarity to train with
lst_CutOffs = [3, 5, 10, 11, 15, 20, 30, 50, 75, 100, 125, 150, 250, 300, 400, 500, 750, 1000, 2000, 5000, 7500, 10000, 15000, 20000, 30000, 33395] # 33395 max. instead of 33397 because of missing descriptors leading to failed KNN

# Loop for all cutoffs
for cutoff in lst_CutOffs:

    neigh = NearestNeighbors(n_neighbors=cutoff)

    # Initialize lists
    ACB_lst_Y_Predictions = []    
    ACB_lst_Y_Original = []
    ACB_lst_amount_training_structures = []
    ACB_lst_max_Similarity = []
    ACB_lst_mean_Similarity = []

    # Loop for all SMILES

    for i in range(0, len(df_Test)):    #

        # Remove test compound from training
        print('Remove SMILES from training:', df_Test['Structure STD'].iloc[i])  
        print('n(Training) BEFORE test cmpnd removal from training data: ', len(df_Train_Pool))
        
        df_Train_Pool_NO_TEST = df_Train_Pool.drop(df_Train_Pool[df_Train_Pool['Structure STD'].isin([df_Test['Structure STD'].iloc[i]])].index, axis='index')


        print('n(Training) AFTER test cmpnd and NA removal from training data: ', len(df_Train_Pool_NO_TEST)) 

        
        print('Calculating distance between test compound and pool...')
        df_Train_Pool_NO_TEST['KNN Distance'] = None

        print('Fitting with NA removal')
        neigh.fit(df_Train_Pool_NO_TEST.loc[:, 'MaxAbsEStateIndex':'fr_urea'].dropna())

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
        ## Find compound cluster with selected tanimoto similarity
        print('Extracting amount of compounds in neighborhood:', cutoff)


        # Initialize lists
        lst_sim_ext_inds = []
        lst_sim_ext_sim_values = []

        try:
            lst_sim_ext_inds = neigh.kneighbors([X_test])[1][0].tolist()           # Index of compounds meeting cutoff (0 - 35k)
            lst_sim_ext_sim_values = neigh.kneighbors([X_test])[0][0].tolist()     # Actual distance values

            print('Highest distance: ', max(lst_sim_ext_sim_values)) 
            print('Lowest distance: ', min(lst_sim_ext_sim_values)) 
        except:
            print('No highest similarity found')  
            lst_sim_ext_sim_values.append(None)   
            lst_sim_ext_inds.append(None)   
            print('Similarities and train indices replaced by "None"')        
    
        
        #_________________________________________________________
        ## DESCRIPTOR CALCULATION + TRAIN

        # Select training and data
        try:           
            print("Extracting 0D RDKit training descriptors from pool for test compound ", i)
            X_train = df_Train_Pool_NO_TEST.loc[:, 'MaxAbsEStateIndex':'fr_urea'].dropna().iloc[lst_sim_ext_inds]   
            y_train = df_Train_Pool_NO_TEST['log10(Papp AB) [cm/s]'].loc[X_train.index]   
        except:
            print('Training descriptor extraction failed for test compound ', i)
            X_train = None
            y_train = None
            print('Train Data replaced by "None"')


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
            ACB_lst_Y_Predictions.append(y_pred[0])  # y_pred (predicted)
            ACB_lst_Y_Original.append(y_test)
        except:    
            print('Appending Pred data failed, filling NA')
            ACB_lst_Y_Predictions.append(None)   # y_pred (predicted)  
            ACB_lst_Y_Original.append(None)  
                                                                             
               
        try:    
            ACB_lst_amount_training_structures.append(len(lst_sim_ext_inds))                # n(Training)
        except:
            print('Appending n(Training) failed, filling NAs')
            ACB_lst_amount_training_structures.append(None) 
        

        try:    
            lst_sim_ext_sim_values_notNone = [x for x in lst_sim_ext_sim_values if x is not None]
            ACB_lst_max_Similarity.append(max(lst_sim_ext_sim_values_notNone))                      # Maximum similarity in training data
        except:    
            print('Appending similarity data failed, filling NAs')            
            ACB_lst_max_Similarity.append(None)                                             # Maximum similarity in training data


        try:    
            # Remove None elements
            lst_sim_ext_sim_values_notNone = [x for x in lst_sim_ext_sim_values if x is not None]
            ACB_lst_mean_Similarity.append(np.mean(lst_sim_ext_sim_values_notNone))                 # Mean similarity in training data
        except:
            ACB_lst_mean_Similarity.append(None)                                            # Mean similarity in training data
       

        print('\n') # New line    

    try:
        # Combine results per cutoff
        ACB_df_Y_test = pd.DataFrame({'Predicted': ACB_lst_Y_Predictions, 
                                    'Measured': ACB_lst_Y_Original,
                                    'n': ACB_lst_amount_training_structures,
                                    'max. Similarity': ACB_lst_max_Similarity,
                                    'mean Similarity': ACB_lst_mean_Similarity})
    except:
        print('Combining results failed.')


    try:    
        # Collect metrics
        var_MSE, var_MAE, var_RMSE, var_R2, var_kTau = calc_regression_metrics(ACB_df_Y_test['Measured'], ACB_df_Y_test['Predicted'])
        dict_metrics = {'MSE': var_MSE, 'MAE': var_MAE, 'RMSE': var_RMSE, 'R2': var_R2, 'K_Tau': var_kTau}
        ACB_df_Metrics = pd.DataFrame([dict_metrics])

        print(ACB_df_Metrics)
        print(ACB_df_Y_test)
    except:
        print('Collecting metrics failed, filling with NAs')
        dict_metrics = {'MSE': None, 'MAE': None, 'RMSE': None, 'R2': None, 'k-tau': None}
        ACB_df_Metrics = pd.DataFrame([dict_metrics])


    try:
        # Export metrics
        file_name_Metrics = 'AR_273x_Metrics_' + str(cutoff) + '_0DRDKit_log10cms_LGBM_altEval_knn'
        ACB_df_Metrics.to_pickle(file_name_Metrics)
    except:
        print('Exporting metrics failed.')


    try:
        # Export results
        file_name_Results = 'AR_273x_Results_' + str(cutoff) + '_0DRDKit_log10cms_LGBM_altEval_knn' # File name dependent on cutoff
        ACB_df_Y_test.to_pickle(file_name_Results) # Save predictions in file
    except: 
        print('Exporting results failed.')



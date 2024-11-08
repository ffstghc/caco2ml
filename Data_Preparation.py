# %% SMILES / MOLECULE STANDARDIZATION AND CLEAN-UP ##
######################################################

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import PandasTools, AllChem, Descriptors, rdMolDescriptors

def standardize_SMILES(dataframe):
    
    lst_inds_Errors_SMILES = []
    smiles_std = []

    for i in range(0, len(dataframe)):
        try:
            standardized = rdMolStandardize.StandardizeSmiles(dataframe['Structure'].iloc[i])
            smiles_std.append(standardized)
        except:
            print("Error SMILES")
            lst_inds_Errors_SMILES.append(i)
            smiles_std.append('Error')
        
    dataframe['Structure STD'] = smiles_std

    return lst_inds_Errors_SMILES

def standardize_MOL(dataframe):
    
    lst_inds_Errors_MOLs = []         
    mol_std = []

    for i in range(0, len(dataframe)):
        try:
            standardized = rdMolStandardize.Cleanup(dataframe['Molecule'].iloc[i])
            mol_std.append(standardized)

        except:
            print("Error MOL")
            lst_inds_Errors_MOLs.append(i)
            mol_std.append('ERROR')  

    dataframe['Molecule STD'] = mol_std 

    return lst_inds_Errors_MOLs

# Standardize SMILES
inds_Errors_SMILES = standardize_SMILES(df) # Add "Structure STD" to dataframe from column "Structure"
df = df.drop(df.index[inds_Errors_SMILES], axis=0) # Drop structures that couldn't be standardized
df = df.dropna()
df = df.reset_index(drop=True)

# Cleanup molecule
PandasTools.AddMoleculeColumnToFrame(df,'Structure STD','Molecule',includeFingerprints=True)
inds_Errors_MOL = standardize_MOL(df) # Add "Molecule STD" to dataframe from column "Molecule"
df = df.drop(df.index[inds_Errors_MOL], axis=0) # Drop molecules that can't be cleaned up
df = df.reset_index(drop=True)
PandasTools.RemoveSaltsFromFrame(df, molCol='Molecule STD') # Remove salts if previous clean-up steps missed it


# %% PADEL DESCRIPTOR CALCULATION ##
####################################

import pandas as pd
from padelpy import from_smiles

def calc_PaDEL_descriptors(dataframe, timeout, n_threads):
    
    dataframe.reset_index(drop=True)
    
    desc = pd.DataFrame()
    
    errors = [] # Initialize list for failed descriptor calculations
    
    for i in range(0, len(dataframe)):
        
        try:
            descriptors = from_smiles(dataframe['Structure STD'].iloc[i], 
                                    fingerprints=False, descriptors=True,
                                    threads=n_threads, timeout=timeout)
        except:
            errors.append(i)     
            
        desc = pd.concat([desc, pd.DataFrame.from_dict([descriptors], orient='columns')])
    
    print("Compound calculated:", i)    
    
    return desc, errors

df_PaDEL_descriptors, lst_indx_PaDEL_errors = calc_PaDEL_descriptors(df, timeout=30, n_threads=16)


# %% CDDD DESCRIPTOR CALCULATION ##
###################################

# CDDD Setup
import pandas as pd
import json
import requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# CDDD REQUEST FUNCTION ##
class CDDDRequest:
    def __init__(self, host="", port=): # Bayer internal API URL & port removed
        self.host = host
        self.port = port
        self.headers = {'content-type': 'application/json'}

    def smiles_to_cddd(self, smiles, preprocess=True):
        url = "{}:{}/smiles_to_cddd/".format(self.host, self.port)
        req = json.dumps({"smiles": smiles, "preprocess": preprocess})
        response = requests.post(url, data=req, headers=self.headers, verify=False)
        return json.loads(response.content.decode("utf-8"))

    def cddd_to_smiles(self, embedding):
        url = "{}:{}/cddd_to_smiles/".format(self.host, self.port)
        req = json.dumps({"cddd": embedding})
        response = requests.post(url, data=req, headers=self.headers, verify=False)
        return json.loads(response.content.decode("utf-8"))

unique_smiles = df['Structure STD']    
cddd_server = CDDDRequest(port=) # Port to Bayer internal CDDD API removed

df = pd.DataFrame(cddd_server.smiles_to_cddd(list(unique_smiles),preprocess=False)) # Convert SMILES to CDDD descriptor



# %% PREPARE CACO-2 DATA #
##########################


## Calcualte log10() column for literature comparison
df_grouped['log10(Papp AB) [cm/s]'] = np.log10(df_grouped['Papp AB [nm/s]'].multiply(10**-7))


# Only include "Reference" assay results
df_grouped.drop(df_grouped[(df_grouped['Result Flag'] != 'Reference') & 
                        (df_grouped['Result Flag'] != 'reference')].index, inplace=True) 


# Only include results from assays with standard 2 umol/L concentration
df_grouped.drop(df_grouped[df_grouped['Conc Test Cmpnd [umol/L]'] != 2].index, inplace=True) 


# Only include results from assays without inhibitor
df_grouped = df_grouped[df_grouped['Conc Inhb 1 [umol/L]'].isnull()]
df_grouped = df_grouped[df_grouped['Conc Inhb 2 [umol/L]'].isnull()]


# Only include compounds with reasonable AB recovery 
df_grouped.drop(df_grouped[ (df_grouped['Recovery AB [%]'] < 50) & 
                        (df_grouped['Recovery AB [%]'] > 200)].index, inplace=True) 

# Only include compounds with reasonable BA recovery 
df_grouped.drop(df_grouped[ (df_grouped['Recovery BA [%]'] < 50) & 
                        (df_grouped['Recovery BA [%]'] > 200)].index, inplace=True) 

# Remove rows with no "Papp AB Passive" information
df_grouped = df_grouped.replace('NA', None).dropna(subset=['Papp AB [nm/s]'])



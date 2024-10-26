# "Exploring the Potential of Adaptive, Local Machine Learning (ML) in Comparison to the Prediction Performance of Global Models: A Case Study from Bayer's Caco-2 Permeability Database"
## _American Chemical Society (ACS): Journal of Chemical Information and Modeling (JCIM)_
### **_Frank Filip Steinbauer, Thorsten Lehr, Andreas Reichel_**
### [DOI PLACEHOLDER]

Repository for archiving the main code chunks used for the local and global machine learning models in the publication **_"Exploring the Potential of Adaptive, Local Machine Learning (ML) in Comparison ton the Prediction Performance of Global Models: A Case Study from Bayer's Caco-2 Permeability Database"_** published in 2024 in **_ACS Journal of Chemical Information and Modeling (JCIM)_** as 1st publication of my doctoral studies at Bayer.

The five different included files contain the main code chunks for:
1. Data preparation (SMILES/molecule object standardization; PaDEL descriptor calculation)
2. Global models (including other descriptor calculations and recursive feature elimination with cross-validation as well as external TDC benchmarking<sup>[1]</sup>)
3. Local model (training data selection via fixed tanimoto similarity criteria)
4. Local model (training data selection via fixed amounts of most similar structuress)
5. Local model (training data selection via kNN as control/proof of superiority of the chosen tanimoto similarity approach)

If you have further questions or need additional parts of the utilized code for your own studies, feel free to contact Filip.Steinbauer@bayer.com.

[1]: https://tdcommons.ai/single_pred_tasks/adme#caco-2-cell-effective-permeability-wang-et-al

#    Random Forest and SHAP for 2D Cellular Metamaterials
#    Johns Hopkins University
#    Shengzhi Luan
#    02.08.2023
#%%  Package Importing

import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import mutual_info_regression

#%%  Subroutine: Mutual Information

def Mutual_Information (feature_select_number, x_norm, y_norm, x):
    
    random_seed = 0
    mi_results = mutual_info_regression(x_norm, y_norm, 
                                        random_state = random_seed)
    header = x.columns.tolist()
    features = header[0:len(header)]
    names_scores = {'Names':features, 'Scores':mi_results}
    ns = pd.DataFrame(names_scores)
    ns = ns.sort_values(by='Scores')
    ns = ns.set_index('Names')
    
    feature_select = ns.sort_values('Scores', ascending = False). \
                                    head(feature_select_number).index.to_list()
    x_select = scaler_x.fit_transform(x[feature_select])
    
    return x_select, feature_select, names_scores

#%% Subroutine: Random Forest Optimal Coefficient

def Coefficient_Searching (estimator_vector, depth_vector, 
                           split_vector, leaf_vector):
    
    random_seed = 42
    model_randF = RandomForestRegressor(random_state=0)
    param_grid = {"n_estimators": estimator_vector, 
                  'max_depth': depth_vector,
                  "min_samples_split": split_vector, 
                  "min_samples_leaf": leaf_vector, 
                  "criterion": ["mse"]}
    
    hp_search = RandomizedSearchCV(estimator = model_randF, 
                                   param_distributions = param_grid, 
                                   n_iter = 100, random_state = random_seed)
    hp_search.fit(x_validation, y_validation)
    
    return hp_search.best_params_

#%% Subroutine: Random Forest Learning

def Random_Forest (n_estimators, max_depth, 
                   min_samples_split, min_samples_leaf, 
                   x_train, y_train, x_validation, y_validation):
    
    random_seed = 42
    model_tree = RandomForestRegressor(random_state = random_seed, 
                                       n_estimators = n_estimators, 
                                       max_depth = max_depth, 
                                       min_samples_split = min_samples_split, 
                                       min_samples_leaf = min_samples_leaf, 
                                       criterion = 'mse')
    model_tree.fit(x_train, y_train)
    y_prediction = model_tree.predict(x_validation)
    model_tree.score(x_validation, y_validation)
    
    return y_prediction, \
           round(model_tree.score(x_validation, y_validation), 3), model_tree

#%% Subroutine: Prediction-Validation Comparison

def Prediction_Plot (y_prediction, y_validation):
    
    data_index = []
    for i in range(len(y_prediction)):
        data_index.append(i)
        
    df_prediction = pd.DataFrame(columns = ['Index', 
                                            'Validation', 
                                            'Prediction'])
    df_prediction['Index'] = data_index
    df_prediction['Validation'] = scaler_y.inverse_transform(y_validation)
    df_prediction['Prediction'] = scaler_y.inverse_transform \
                                  (y_prediction.reshape(len(y_prediction),-1))
    
    fig = plt.figure(figsize = (15,10))
    plt.scatter(df_prediction['Validation'], df_prediction['Prediction'])
    plt.plot([0,0.1],[0,0.1],'k--')
    plt.plot([], [], ' ', label = "R-squared: " + str(score))
    plt.legend()
    plt.xlabel('Validation', fontsize = 18)
    plt.ylabel('Prediction', fontsize = 18)
    
    return fig, df_prediction

#%% Subroutine: SHAP General Learning

def SHAP_Plot (shap_values, samples, feature_select):
    
    shap.summary_plot(shap_values, samples, feature_names = feature_select, 
                      show = False, plot_type = 'violin', plot_size = (10,10), 
                      max_display = 20)
    plt.savefig('SHAP_Figure.png', dpi = 500)
    
#%% Subroutine: SHAP Individual Learning

def SHAP_Force_Plot (sample_feature, sample_value, sample_index):
    
    all_feature = sample_feature.append(x)
    all_feature_select = scaler_x.fit_transform(all_feature[feature_select])
    sample_feature_select = all_feature_select[0:len(sample_feature)]
    sample_prediction_norm = model_tree.predict(sample_feature_select)
    sample_prediction = scaler_y.inverse_transform \
                        (sample_prediction_norm.reshape \
                         (len(sample_prediction_norm), -1))
    individual_result = pd.DataFrame(columns = ['Validation', 'Prediction'])
    individual_result['Validation'] = sample_value
    individual_result['Prediction'] = sample_prediction
    
    sample_SHAP = explainer.shap_values(sample_feature_select, 
                                        approximate = False, 
                                        check_additivity = False)
    shap.force_plot(explainer.expected_value, sample_SHAP[sample_index-1,:], 
                    np.around(sample_feature_select[sample_index-1,:], 
                              decimals = 2), 
                    matplotlib = True, show = False, 
                    feature_names = feature_select, figsize = (20,3), 
                    contribution_threshold = 0.15)
    
    return individual_result

#%% Importing Data

data_path = 'D:/Johns Hopkins University/Research/2023/January-March/Machine Learning/Random Forest/Structure-Property-Refined-Data.csv'
data_matrix = pd.read_csv(data_path)
x = data_matrix[data_matrix.columns[0:len(data_matrix.columns)-1]]
y = data_matrix[data_matrix.columns[len(data_matrix.columns)-1]]

#%% Hyperparameter

feature_select_number = 23
estimator_vector = [100,200,300,400,500]
depth_vector = [10,20,30,40,50]
split_vector = [0.1,0.01,0.001,0.0001,0.00001]
leaf_vector = [1,2,3,4,5]

#%% Data Normalization

scaler_x = MinMaxScaler(feature_range = (0,1))
scaler_y = MinMaxScaler(feature_range = (0,1))
x_norm = scaler_x.fit_transform(x)
y_norm = scaler_y.fit_transform((y.to_numpy()).reshape(-1,1))

#%% Random Forest Machine Learning

x_select, feature_select, names_scores = Mutual_Information (feature_select_number, x_norm, y_norm, x)
x_select = pd.DataFrame(x_select)
y_norm = pd.DataFrame(y_norm)
x_train, x_validation , y_train, y_validation = train_test_split(x_select, y_norm, test_size = 0.30, random_state = 79)
# rf_coefficient = Coefficient_Searching(estimator_vector, depth_vector, split_vector, leaf_vector)
# input_coefficient = [rf_coefficient['n_estimators'], rf_coefficient['max_depth'], rf_coefficient['min_samples_split'], rf_coefficient['min_samples_leaf']]
input_coefficient = [1000, 50, 0.001, 5]
y_prediction, score, model_tree = Random_Forest(input_coefficient[0], input_coefficient[1], input_coefficient[2], input_coefficient[3], x_train, y_train, x_validation, y_validation)
fig, df_prediction = Prediction_Plot(y_prediction, y_validation)

#%% SHAP General Analysis

shap.initjs()
explainer = shap.TreeExplainer(model_tree)
shap_values = explainer(x_validation)
shap_values = explainer.shap_values(x_validation, approximate = False, check_additivity = False)
SHAP_Plot(shap_values, x_validation, feature_select)

#%% SHAP Individual Analysis

index_of_sample = [0,50]
sample_index = 1
sample_matrix = data_matrix.iloc[index_of_sample]
sample_feature = sample_matrix[sample_matrix.columns[0:len(sample_matrix.columns)-1]]
sample_value = sample_matrix[sample_matrix.columns[len(sample_matrix.columns)-1]]
SHAP_Force_Plot (sample_feature, sample_value, sample_index)
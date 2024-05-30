# %%
from sklearn.datasets import make_regression
from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
import pandas as pd

# Generate synthetic data
X, y = make_regression(n_samples=100, n_features=20, n_targets=3, random_state=1)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

# Initialize the StandardScaler
scaler = StandardScaler()

# Fit the scaler on the training data and transform both training and testing data
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Initialize the StandardScaler for targets if you wish to scale them
target_scaler = StandardScaler()
# Fit the scaler on the training target data and transform both training and testing target data
y_train_scaled = target_scaler.fit_transform(y_train)
y_test_scaled = target_scaler.transform(y_test)
# Initialize the SVM regressor
svm_regressor = SVR(kernel='linear')
# Initialize the Ridge regressor
ridge_regressor = Ridge(alpha=1.0)

# Create the MultiOutputRegressor
multi_target_regressor = MultiOutputRegressor(ridge_regressor)

# Fit the model on the training data
multi_target_regressor.fit(X_train, y_train)

# Make predictions on the test data
predictions = multi_target_regressor.predict(X_test)

# Example usage:
# Here we take the first test example and predict its targets
example = X_test[0]
predicted_targets = multi_target_regressor.predict([example])
print(f"Predicted targets: {predicted_targets[0]} actual {y_test[0]}")

original_scale_predictions = target_scaler.inverse_transform(predictions)

# %%
df = pd.read_csv(r"H:\GitRepos\ShipYard2\data\PMM-results\PMM_results .csv",header=0,delimiter=';')
#remove unnamed columns
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
df.dropna(how='any', inplace=True)  # Drops rows with any NaN value

# create the input and output 
# %%
X = df[['Lpp','Tmean/Lpp','B/Tmean','LCG/Lpp','Fn','Pdiam/Tmean','RudArea/(Lpp*Tmean)','Trim=(Taft-Tfore)/Lpp','BlockCoefficient']]
y_xCoefs=df[['xv','xvv','xrr','xvr','xudot']]
y_yCoefs=df[['yv','yvv','yvabsv','yr','yru','yrrr','yvrr','yrvv','yrabsv','yvdot','yrdot','yd']]
y_nCoefs=df[['nv','nvabsv','nr','nrrr','nvrr','nrvv','nrabsv','nvabsr','nvdot','nrdot','nvv','nvu','nru','nd']]

X_test=X.loc[[4,42]]
X_train=X.drop([4,42])

y_xtest=y_xCoefs.loc[[4,42]]
y_xtrain=y_xCoefs.drop([4,42])

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

y_ytest=y_yCoefs.loc[[4,42]]
y_ytrain=y_yCoefs.drop([4,42])

y_ntest=y_nCoefs.loc[[4,42]]
y_ntrain=y_nCoefs.drop([4,42])
# %%
##########
# Fit the model on the training data
svm_regressor = SVR(kernel='linear')
# Initialize the Ridge regressor
ridge_regressor = Ridge(alpha=1.0)

# Create the MultiOutputRegressor
multi_target_regressor = MultiOutputRegressor(ridge_regressor)

multi_target_regressor.fit(X_train, y_xtrain)
predictions = multi_target_regressor.predict(X_test)
example = X_test.loc[4]
predicted_targets = multi_target_regressor.predict([example])
print(f"Predicted :      actual ")
print(f'xv {predicted_targets[0][0]:.1f}        {y_xtest.iloc[0,0]}')
print(f'xvv {predicted_targets[0][1]:.1f}       {y_xtest.iloc[0,1]}')
print(f'xrr {predicted_targets[0][2]:.1f}       {y_xtest.iloc[0,2]}')
print(f'xvr {predicted_targets[0][3]:.1f}       {y_xtest.iloc[0,3]}')
print(f'xudot {predicted_targets[0][4]:.1f}     {y_xtest.iloc[0,4]}')

# %%

# %%
from sklearn.datasets import make_regression
from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split

# Generate synthetic data
X, y = make_regression(n_samples=100, n_features=20, n_targets=3, random_state=1)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

# Initialize the base regressor
ridge = Ridge()

# Create the MultiOutputRegressor
multi_target_regressor = MultiOutputRegressor(ridge)

# Fit the model on the training data
multi_target_regressor.fit(X_train, y_train)

# Make predictions on the test data
predictions = multi_target_regressor.predict(X_test)

# Example usage:
# Here we take the first test example and predict its targets
example = X_test[0]
predicted_targets = multi_target_regressor.predict([example])
print(f"Predicted targets: {predicted_targets[0]}")

# %%

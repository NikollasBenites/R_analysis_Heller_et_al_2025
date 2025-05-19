library(ggplot2)
library(caret)


# Linear Regression R base ------------------------------------------------

plot_lr_base <- function(data, x, y, 
                                   title = "Linear Regression", 
                                   xlab = NULL, ylab = NULL) {
  # Convert to dataframe
  data <- as.data.frame(data)
  X <- data[[x]]
  Y <- data[[y]]
  
  # Fit the linear model
  formula <- as.formula(paste("`", y, "` ~ `", x, "`", sep = ""))
  lr <- lm(formula, data = data)
  r_sqrd <- summary(lr)$r.squared
  
  # Set axis labels (default to column names if NULL)
  x_label <- if (!is.null(xlab)) xlab else x
  y_label <- if (!is.null(ylab)) ylab else y
  
  # Create scatter plot with regression line
  plot(X, Y, main = title, xlab = x_label, ylab = y_label, pch = 16, col = "black")
  abline(lr, col = "red", lwd = 2)
  
  # Add R² to the top-right corner
  text(max(X), max(Y), labels = paste("R² = ", round(r_sqrd, 2)), 
       pos = 2, col = "black", font = 2, cex = 1, adj = c(1, 1))
  
  # Print model summary
  print(summary(lr))
}

# -------------------- BASE LINEAR REGRESSION FUNCTION LAYOUT --------------------
base_lr_lo <- function(data, x_col, y_col) {
  # Convert to dataframe
  data <- as.data.frame(data)
  X <- data[[x_col]]
  Y <- data[[y_col]]
  
  # Fit a linear model with corrected formula syntax
  formula <- as.formula(paste("`", y_col, "` ~ `", x_col, "`", sep = ""))
  lr <- lm(formula, data = data)
  r_sqrd <- summary(lr)$r.squared
  
  # Create a 2x2 plot layout
  par(mfrow = c(2, 2))
  
  # Scatter plot with regression line
  plot(X, Y, main = "Linear Regression (Base)", xlab = x_col, ylab = y_col, pch = 1)
  abline(lr, col = "red", lwd = 1)
  
  # Add R² to top-right corner
  text(max(X), max(Y), labels = paste("R² = ", round(r_sqrd, 2)), 
       pos = 2, col = "black", font = 2, cex = 0.9, adj = c(1, 1))
  
  # Diagnostic plots
  plot(lr, which = 1)  # Residuals vs Fitted
  plot(lr, which = 2)  # Normal Q-Q
  plot(lr, which = 3)  # Scale-Location
  
  # Print model summary
  print(summary(lr))
}

# -------------------- CARET LINEAR REGRESSION FUNCTION LAYOUT--------------------
caret_lr_lo <- function(data, x_col, y_col) {
  # Convert to dataframe
  data <- as.data.frame(data)
  
  set.seed(123)  # Ensure reproducibility
  train_control <- trainControl(method = "cv", number = 10)
  
  # Train the linear model using caret (with correct column name handling)
  formula <- as.formula(paste("`", y_col, "` ~ `", x_col, "`", sep = ""))
  lr_model <- train(formula, data = data, method = "lm", trControl = train_control)
  
  # Extract final model coefficients & R²
  coef <- coef(lr_model$finalModel)
  intercept <- coef[1]
  slope <- coef[2]
  r_sqrd <- summary(lr_model$finalModel)$r.squared
  
  # Create a 2x2 plot layout
  par(mfrow = c(2, 2))
  
  # Scatter plot with regression line
  plot(data[[x_col]], data[[y_col]], main = "Linear Regression (Caret)", 
       xlab = x_col, ylab = y_col, pch = 1)
  abline(intercept, slope, col = "red", lwd = 1)
  
  # Add R² to top-right corner
  text(max(data[[x_col]]), max(data[[y_col]]), labels = paste("R² = ", round(r_sqrd, 2)), 
       pos = 2, col = "black", font = 2, cex = 0.9, adj = c(1, 1))
  
  # Diagnostic plots
  plot(lr_model$finalModel, which = 1)  # Residuals vs Fitted
  plot(lr_model$finalModel, which = 2)  # Normal Q-Q
  plot(lr_model$finalModel, which = 3)  # Scale-Location
  
  # Print model summary
  print(lr_model)
  print(summary(lr_model$finalModel))
}

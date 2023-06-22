# Change the working directory
setwd('/Users/luisquinonespr/code/BSE/mlfinance/hw2')

# Load the RDS file
data <- readRDS('SP500.rds')

# Write to a CSV file
write.csv(data, 'SP500.csv')

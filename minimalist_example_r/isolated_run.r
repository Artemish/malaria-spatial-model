source("train.r")
source("predict.r")

train_chap(
  "input/training_data_harmonized.csv",
  "output/model.bin"
)

predict_chap(
  "output/model.bin",
  "input/training_data_harmonized.py",
  "input/test_data_harmonized.csv",
  "output/predictions.csv"
)

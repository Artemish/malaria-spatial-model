import sys
import pandas as pd

formatted_outfile='training_data_formatted.csv'
rscript = ""

def reshape_input(df):

if __name__ == '__main__':
    data_input_path = sys.argv[1]
    input_data = pd.read_csv(data_input_path)

    input_data_reformatted = reshape_input(input_data)
    input_data_reformatted.to_csv(formatted_outfile)

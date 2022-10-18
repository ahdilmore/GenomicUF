import pandas as pd
import glob
import os
METADATA_COLUMNS = ['mouse_name', 'time_point', 'isolate_number',
                    'sequencing_lane']

def _split_filename(x, col_name):
    name_list = x.split('_')
    if col_name == 'mouse_name':
        return name_list[0]
    elif col_name == 'time_point':
        return name_list[1]
    elif col_name == 'isolate_number':
        return name_list[2]
    elif col_name == 'sequencing_lane':
        return name_list[3]

def make_metadata(sample_names):
    df = pd.DataFrame(data={'sample_name': sample_names})
    for colname in METADATA_COLUMNS:
        df.insert(loc=0, column=colname,
                  value=df['sample_name'].apply(_split_filename,
                  col_name=colname))
    df = df.set_index('sample_name')
    return df

def main():
    bacteria_dict = {'AZ20': 'adaptation_AZ20/', 'AZ51': 'adaptation_AZ51/'}

    for bacteria in bacteria_dict:
        if not os.path.exists(bacteria_dict[bacteria] + 'metadata.csv'): 
            samples = [os.path.basename(x) for x in glob.glob(bacteria_dict[bacteria] + '*')]
            md_df = make_metadata(samples)
            md_df.insert(loc=0, column='bacteria', value=bacteria)
            md_df.to_csv(bacteria_dict[bacteria] + 'metadata.csv')

if __name__ == "__main__":
    main()
import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
##################################################################
def MatchingFunction(pattern, string, number_of_match):
        #'''return array with the matchs'''
        output_array = []
        match = re.search(pattern,string)
        if match:
            i = 0
            while i < int(number_of_match):
                temp = match.group(i+1)
                output_array.append(temp)
                i= i + 1
            return output_array 
        elif not match:
            #print ('Nothing match - ', string)
            return 0
##################################################################
def worker(string, pattern, file_name):
    number_of_match = 6
    matched_array =  MatchingFunction(pattern, string, number_of_match)
    #print(matched_array[0], matched_array[1], matched_array[2], matched_array[3], matched_array[4], matched_array[5])
    return  matched_array[0], matched_array[1], matched_array[2], matched_array[3], matched_array[4], matched_array[5], file_name
##################################################################
def EnhancerNameCreator(chrom, start, end):
    string = chrom+'_'+str(start)+'_'+str(end)
    return string
##################################################################
if __name__ == '__main__':

    folder1 = 'ChiA_PET'
    folder2 = 'pol2_chia_pet_Mcf7'
    file_name_array = []
    file1 = 'CHM053T_hg19.ChromatinInteractons.bed'
    file2 = 'CHM040M_hg19.ChromatinInteractons.bed'
    file3 = 'CHM160M_L4.ChromatinInteractons.bed'
    file4 = 'CHM163M_L4.ChromatinInteractons.bed'
    file_name_array.append(file1)
    file_name_array.append(file2)
    file_name_array.append(file3)
    file_name_array.append(file4)

    paths_array = []
    for file in file_name_array:
        path = os.path.join(folder1,folder2,file)
        print(path)
        paths_array.append(path)

    df_array = []
    for path in paths_array:
        df1 = pd.read_csv(path, sep = '\t', header = None)
        print(len(df1))
        #print(df1.head(3))
        df_array.append(df1)
    print(len(df_array))

    df_new_array = []
    for df_source in df_array:
        df = pd.DataFrame(index = range(len(df_source)), columns = ['CHR1','Start1','End1','CHR2','Start2','End2','File_Name'])
        df_new_array.append(df)
    print(len(df_new_array))

    full_DF_ChiA_Pet = pd.DataFrame()
    pattern = re.compile('chr(\w+)\:(\d+)\.\.(\d+)\-chr(\w+)\:(\d+)\.\.(\d+)\,\d+',re.M|re.I)
    for df, df_source, file_name in zip(df_new_array, df_array, file_name_array):
        df['CHR1'],df['Start1'],df['End1'], df['CHR2'],df['Start2'],df['End2'], df['File_Name'] =\
        np.vectorize(worker)(df_source[3], pattern, file_name)
        print(len(df))
        full_DF_ChiA_Pet = pd.concat([full_DF_ChiA_Pet, df])

    full_DF_ChiA_Pet.drop_duplicates(inplace = True)
    full_DF_ChiA_Pet.reset_index(drop=True, inplace = True)
    full_DF_ChiA_Pet['Code'] = full_DF_ChiA_Pet.index

    df1_chia_pet = full_DF_ChiA_Pet[['CHR1','Start1', 'End1','File_Name', 'Code']]
    df2_chia_pet = full_DF_ChiA_Pet[['CHR2','Start2', 'End2','File_Name','Code']]

    if not os.path.exists('output'):
        os.makedirs('output')

    df1_chia_pet.to_csv('output/ChIA_PET_Side1.csv',sep = ',', index = False)
    df2_chia_pet.to_csv('output/ChIA_PET_Side2.csv',sep = ',', index = False)









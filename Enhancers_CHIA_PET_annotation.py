import pandas as pd
import numpy as np
import re
import os
import glob
import multiprocessing

######################################################################
def EnhancerInput():
    folderup = '..'
    folder1 = 'Step2_Rebin_Enhancers'
    folder2 = 'Enhancers_list'
    filename = 'Enhancers_list_reduced.csv'
    path = os.path.join(folderup,folder1,folder2,filename)
    df_en = pd.read_csv(path, sep = ',',low_memory=False)
    return df_en
######################################################################
def ConductAreasInput():
    folderup = '..'
    folder1 = 'Step8_ChIA_PET_Sides_Separation'
    folder2 = 'output'
    filename1 = 'ChIA_PET_Side1.csv'
    filename2 = 'ChIA_PET_Side2.csv'
    path1 = os.path.join(folderup,folder1,folder2,filename1)
    path2 = os.path.join(folderup,folder1,folder2,filename2)
    df1 = pd.read_csv(path1, sep = ',', low_memory = False)
    df2 = pd.read_csv(path2, sep = ',', low_memory = False)
    return df1, df2
######################################################################
def MakeDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
######################################################################
def SaveDFtoPickle(directory, output_name, df):
    MakeDirectory(directory)
    path = os.path.join(directory, str(output_name))
    df.to_pickle(path)
        
######################################################################
def SaveDFtoCSV(directory,output_name, df):
    MakeDirectory(directory)
    path = os.path.join(directory, str(output_name))
    df.to_csv(path, sep=',', index = False)
######################################################################
def MultiprocessFunctionTreeInput(function, df_chrom_array, df1_array, df2_array):
    results = []
    #results_async = [MultiprocessingChrom(chrom, df_1_c, df_2_c)\
    # for chrom, df_1_c, df_2_c in zip(df_chrom_array, df1_array, df2_array)]
    results_async = [pool.apply_async(function, (chrom, df_1_c, df_2_c))\
     for chrom, df_1_c, df_2_c in zip(df_chrom_array, df1_array, df2_array)]
    results=[r.get() for r in results_async]

    return results
######################################################################
def FinalDataFrameReconstruction(df_array):
    df_final = pd.DataFrame()
    for df1 in df_array:
        df_final = pd.concat([df_final, df1])
    print('Final DataFrame:\n{0}\n'.format(df_final.head()))
    return df_final
######################################################################
def LoopingWorker(df_en, df_chia):
    df_chia.columns = ['CHR', 'Start','End', 'File_Name', 'Code']
    chrom_list = set(df_chia['CHR'])
    #print(chrom_list)
    #print(chrom_list1)
    df_chia_array = []
    df_en_array = []
    df_chrom_array = []
    division_n = 0

    for chrom in chrom_list:
        df_chia_chrom = df_chia[df_chia['CHR'] == str(chrom)]
        chrom_en = 'chr'+ str(chrom)
        df_en_chrom = df_en[df_en['CHR'] == chrom_en]
        print('The number of position at the chrom {0} is {1}'.format(chrom, len(df_chia_chrom)))
        print('The number of Enhancers at the chrom {0} is {1}'.format(chrom, len(df_en_chrom)))
        
        division_n = (len(df_chia_chrom)//500)
        if division_n == 0:
            df_chia_array.append(df_chia_chrom)
            df_en_array.append(df_en_chrom)
            df_chrom_array.append(chrom)
            #print(division_n, len(df_450_chrom), len(df_en_chrom), chrom)
        elif division_n > 0:
            df1_temp_list = np.array_split(df_chia_chrom, division_n)
            for df in df1_temp_list:
                df_chia_array.append(df)
                df_en_array.append(df_en_chrom)
                df_chrom_array.append(chrom)
                #print(division_n, len(df), len(df_en_chrom), chrom)
        else:
            print('Why here?')

        #break
        #print (chrom, len(df_450_array), len(df_en_array), len(df_chrom_array ))
    return df_chrom_array, df_chia_array, df_en_array
######################################################################
def DesisionMaking(start_en, end_en, start_chia, end_chia):
    if int(start_en) <= int(start_chia) and int(start_chia) <= int(end_en):
        return 1
    elif int(end_chia) >= int(start_en) and int(end_chia) <= int(end_en):
        return 1
    elif int(start_chia) <= int(start_en) and int(end_chia) >= int(end_en):
        return 1
    else:
        return 0
######################################################################
def MultiprocessingChrom(chrom, df_chia_c, df_en_c):
    ######################################################################
    def SecondRound(start_chia, end_chia, df_en_c = df_en_c):
        df_en_c['Pos'] = np.vectorize(DesisionMaking)(df_en_c['Start'], df_en_c['End'], start_chia, end_chia)
        df_en_c_pos = df_en_c[df_en_c['Pos'] == 1]
        if len(df_en_c_pos) > 1:
            array_t = list(df_en_c_pos['Enhancer_Name'])
            myString = ",".join(array_t)
            return myString
        elif len(df_en_c_pos) == 1:
            array_t = list(df_en_c_pos['Enhancer_Name'])
            return array_t[0]
        else:
            return 'None'
    ######################################################################
    #df_chia_c.columns = ['CHR','Start','End', 'File_Name', 'Code']
   
    df_chia_c['Annotation_En'] = np.vectorize(SecondRound)(df_chia_c['Start'], df_chia_c['End'])

    return df_chia_c
    
######################################################################
def EnhancerName(chrom, start, end):
    mystring = chrom+'_'+str(start)+'_'+str(end)
    return mystring

######################################################################
if __name__ == '__main__':
    pool = multiprocessing.Pool()
    df_en  = EnhancerInput()
    df_en['Enhancer_Name'] = np.vectorize(EnhancerName)(df_en['CHR'], df_en['Start'], df_en['End'])
    print(len(df_en))
    print(df_en.head())
    df_chia1, df_chia2 = ConductAreasInput()
    print(len(df_chia1))
    print(len(df_chia2))
    print(df_chia1.head())
    print(df_chia2.head())
    df_chrom_array1, df_chia_array1, df_en_array1 = LoopingWorker(df_en, df_chia1)
    results1 = MultiprocessFunctionTreeInput(MultiprocessingChrom, df_chrom_array1, df_chia_array1, df_en_array1)
    final_df1 = FinalDataFrameReconstruction(results1)
    SaveDFtoCSV('Results','Annotated_Chia_List1.csv', final_df1)

    df_chrom_array2, df_chia_array2, df_en_array2 = LoopingWorker(df_en, df_chia2)
    results2 = MultiprocessFunctionTreeInput(MultiprocessingChrom, df_chrom_array2, df_chia_array2, df_en_array2)
    final_df2 = FinalDataFrameReconstruction(results2)
    SaveDFtoCSV('Results','Annotated_Chia_List2.csv', final_df2)




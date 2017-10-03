import pandas as pd
import numpy as np
import re
import os
import glob
import multiprocessing

######################################################################
def GenesInput():
    folderup = '..'
    folder1 = 'Step1_Remove_Enhancer_Near_TSS'
    filename = 'gh19_RefSeq_Genes.txt'
    path = os.path.join(folderup,folder1,filename)
    df = pd.read_csv(path,  sep ='\t', error_bad_lines=False)
    df_gene_small = df[['name','chrom','strand','txStart','txEnd','name2']]
    return df_gene_small
######################################################################
def ConductAreasInputAnnotatedwithEnhancers():
    folderup = '..'
    folder1 = 'Step9_Enhancers_ChIA_PET_annotation'
    folder2 = 'Results_Supper_Enhancers'
    filename1 = 'Annotated_Chia_List1.csv'
    filename2 = 'Annotated_Chia_List2.csv'
    path1 = os.path.join(folderup,folder1,folder2,filename1)
    path2 = os.path.join(folderup,folder1,folder2,filename2)
    df1 = pd.read_csv(path1, sep = ',', low_memory = False)
    df2 = pd.read_csv(path2, sep = ',', low_memory = False)
    return df1, df2
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
def LoopingWorker(df_g, df_chia):
    df_chia.columns = ['CHR', 'Start','End', 'File_Name', 'Code']
    chrom_list = list(range(1,22)) + ['X','Y']
    print(chrom_list)
    #print(chrom_list1)
    df_chia_array = []
    df_g_array = []
    df_chrom_array = []
    division_n = 0

    for chrom in chrom_list:
        print(chrom)
        df_chia_chrom = df_chia[df_chia['CHR'] == str(chrom)]
        chrom_en = 'chr'+ str(chrom)
        df_g_chrom = df_g[df_g['CHR'] == chrom_en]
        print('The number of position at the chrom {0} is {1}'.format(chrom, len(df_chia_chrom)))
        print('The number of Genes at the chrom {0} is {1}'.format(chrom, len(df_g_chrom)))

        division_n = (len(df_chia_chrom)//500)
        if division_n == 0:
            df_chia_array.append(df_chia_chrom)
            df_g_array.append(df_g_chrom)
            df_chrom_array.append(chrom)
            #print(division_n, len(df_450_chrom), len(df_en_chrom), chrom)
        elif division_n > 0:
            df1_temp_list = np.array_split(df_chia_chrom, division_n)
            for df in df1_temp_list:
                df_chia_array.append(df)
                df_g_array.append(df_g_chrom)
                df_chrom_array.append(chrom)
                #print(division_n, len(df), len(df_en_chrom), chrom)
        else:
            print('Why here?')

        #break
        #print (chrom, len(df_450_array), len(df_en_array), len(df_chrom_array ))
    return df_chrom_array, df_chia_array, df_g_array
######################################################################
def DesisionMaking(TSS, start_chia, end_chia):
    dif = 1500
    start_chia = int(start_chia) - dif
    end_chia = int(end_chia) + dif

    if int(start_chia) <= int(TSS) and int(TSS) <= int(end_chia):
        return 1
    else:
        return 0
######################################################################
def MultiprocessingChrom(chrom, df_chia_c, df_g_c):

    ######################################################################
    def SecondRound(start_chia, end_chia, df_en_c = df_g_c):
        if len(df_en_c) > 0:
            df_g_c['Pos'] = np.vectorize(DesisionMaking)(df_g_c['TSS'], start_chia, end_chia)
            df_g_c_pos = df_g_c[df_g_c['Pos'] == 1]
            if len(df_g_c_pos) > 1:
                array_t = list(df_g_c_pos['Gene_ID'])
                array_t1 = list(df_g_c_pos['Gene_Name'])
                array_t2 = list(df_g_c_pos['TSS'])
                myString = ",".join(array_t)
                myString1 = ",".join(array_t1)
                myString2 = ','.join(str(x) for x in array_t2)
                return myString, myString1, myString2
            elif len(df_g_c_pos) == 1:
                array_t = list(df_g_c_pos['Gene_ID'])
                array_t1 = list(df_g_c_pos['Gene_Name'])
                array_t2 = list(df_g_c_pos['TSS'])
                return array_t[0], array_t1[0], str(array_t2[0])
            else:
                return 'None', 'None', 'None'
        else:
            return 'None', 'None', 'None'
    ######################################################################
    #df_chia_c.columns = ['CHR','Start','End', 'File_Name', 'Code']
    if len(df_chia_c) > 0:
        df_chia_c['Gene_ID'], df_chia_c['Gene_Name'], df_chia_c['TSS']  = np.vectorize(SecondRound)(df_chia_c['Start'], df_chia_c['End'])

    return df_chia_c
    
######################################################################
def NegativeDNAStrandCorrection(gene_strand, gene_start, gene_end):
    if gene_strand == '+':
        return gene_start
    elif gene_strand == '-':
        return gene_end
    else:
        print('Error on TSSs')
######################################################################
if __name__ == '__main__':
    pool = multiprocessing.Pool()
    
    df_genes  = GenesInput()
    df_genes.columns = ['Gene_ID', 'CHR', 'strand', 'txStart', 'txEnd', 'Gene_Name']
    df_genes['TSS'] = np.vectorize(NegativeDNAStrandCorrection)(df_genes['strand'], df_genes['txStart'], df_genes['txEnd'])
    df_g = df_genes[['Gene_ID', 'Gene_Name', 'CHR', 'TSS']]
    print(len(df_g))
    print(df_g.head())
    
    df_chia1, df_chia2 = ConductAreasInput()
    print(len(df_chia1))
    print(len(df_chia2))
    print(df_chia1.head())
    print(df_chia2.head())

    df_chia1_en, df_chia2_en = ConductAreasInputAnnotatedwithEnhancers()
    print(df_chia1_en.head())
    print(df_chia2_en.head())

    df_chrom_array1, df_chia_array1, df_g_array1 = LoopingWorker(df_g, df_chia1)
    results1 = MultiprocessFunctionTreeInput(MultiprocessingChrom, df_chrom_array1, df_chia_array1, df_g_array1)
    final_df1 = FinalDataFrameReconstruction(results1)
    final_df1_1 = pd.merge(df_chia1_en, final_df1, on = ['CHR', 'Start', 'End', 'File_Name', 'Code'], how = 'outer')
    SaveDFtoCSV('Results_Super_Enhancers','Gene_Annotated_Chia_List1.csv', final_df1_1)

    df_chrom_array2, df_chia_array2, df_g_array2 = LoopingWorker(df_g, df_chia2)
    results2 = MultiprocessFunctionTreeInput(MultiprocessingChrom, df_chrom_array2, df_chia_array2, df_g_array2)
    final_df2 = FinalDataFrameReconstruction(results2)
    final_df2_2 = pd.merge(df_chia1_en, final_df1, on = ['CHR', 'Start', 'End', 'File_Name', 'Code'], how = 'outer')
    SaveDFtoCSV('Results_Super_Enhancers','Gene_Annotated_Chia_List2.csv', final_df2_2)




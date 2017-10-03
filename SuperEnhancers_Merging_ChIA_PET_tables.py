import pandas as pd
import numpy as np
import os

######################################################################
def MakeDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
######################################################################
def SaveDFtoCSV(directory,output_name, df):
    MakeDirectory(directory)
    path = os.path.join(directory, str(output_name))
    df.to_csv(path, sep=',', index = False)
##############################################################
def PairingAnnotation(annot_en1, annot_gene1, annot_en2, annot_gene2):
    enhancer_promoter = 'None'
    enhancer_enhancer = 'None'
    promoter_promoter = 'None'
    if annot_en1 != 'None' and annot_gene2 != 'None':
        enhancer_promoter = 'Yes'
    if annot_en2 != 'None' and annot_gene1 != 'None':
        enhancer_promoter = 'Yes'
    if annot_en1 != 'None' and annot_en1 != 'None':
        enhancer_enhancer = 'Yes'
    if annot_gene2 != 'None' and annot_gene1 != 'None':
        promoter_promoter = 'Yes'

    return enhancer_promoter, enhancer_enhancer, promoter_promoter
##############################################################
def ReOrdering(Annotation_En_1,Gene_ID_1,Gene_Name_1,Annotation_En_2,Gene_ID_2,Gene_Name_2):
    if Annotation_En_1 == 'None' and Annotation_En_2 != 'None':
        return Annotation_En_2, Gene_ID_2,Gene_Name_2, Annotation_En_1,Gene_ID_1,Gene_Name_1
    else:
        return Annotation_En_1,Gene_ID_1,Gene_Name_1,Annotation_En_2,Gene_ID_2,Gene_Name_2

##############################################################
def DeconvolveDataFrame(df, colname_split, colname_code):
    df = df.fillna('None')
    
    df1 = pd.DataFrame(df[colname_split].str.split(',').tolist(), index = df[colname_code]).stack()
    df1 = df1.reset_index()[[0, colname_code]]
    df1.columns = [colname_split, colname_code]
    
    df_n = df.drop([colname_split], axis =1)
    
    df_m = pd.merge(df_n, df1, on = colname_code, how = 'outer')
    
    return df_m
##############################################################
def ReorderingPromoter(Annotation_En_1, Gene_ID_1, Annotation_En_2, Gene_ID_2):
    
    if Gene_ID_1 != 'None' and Gene_ID_2 == 'None':
        Gene_ID_2 = Gene_ID_1
        Gene_ID_1 = 'None'

    if Annotation_En_1 == 'None'  and Annotation_En_2 != 'None':
        Annotation_En_1 = Annotation_En_2
        Annotation_En_2 = 'None'

    return Annotation_En_1, Gene_ID_1, Annotation_En_2, Gene_ID_2
##############################################################
if __name__ == '__main__':
    folder_up = '..'
    folder1 = 'Step10_Gene_ChIA_PET_annotation'
    folder2 = 'Results_Super_Enhancers'
    filename1 = 'Gene_Annotated_Chia_List1.csv'
    filename2 = 'Gene_Annotated_Chia_List2.csv'
    path1 = os.path.join(folder_up, folder1,folder2, filename1)
    path2 = os.path.join(folder_up, folder1, folder2, filename2)

    df1 = pd.read_csv(path1, sep = ',')
    df2 = pd.read_csv(path2, sep = ',')

    print(df1.head())
    print(df1.head())
    
    df_final = pd.merge(df1,df2, on = 'Code', how = 'outer', suffixes = ('_1','_2'))
    
    df_final['En-Pro'],df_final['En-En'] ,df_final['Pro-Pro']\
     = np.vectorize(PairingAnnotation)(df_final['Annotation_En_1'], df_final['Gene_ID_1'],df_final['Annotation_En_2'], df_final['Gene_ID_2'])

    df_en_pro = df_final[df_final['En-Pro'] == 'Yes']
    df_pro_pro = df_final[df_final['Pro-Pro'] == 'Yes']
    df_en_en = df_final[df_final['En-En'] == 'Yes']

    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated.csv',df_final)
    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_en_pro.csv',df_en_pro)
    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_pro_pro.csv',df_pro_pro)
    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_en_en.csv',df_en_en)

    df_en_pro_small = df_en_pro[['Annotation_En_1','Gene_ID_1','Gene_Name_1','Annotation_En_2','Gene_ID_2','Gene_Name_2']]
    df_en_pro_rearranged = pd.DataFrame()
    df_en_pro_rearranged['Annotation_En_1'], df_en_pro_rearranged['Gene_ID_1'], df_en_pro_rearranged['Gene_Name_1'],\
    df_en_pro_rearranged['Annotation_En_2'], df_en_pro_rearranged['Gene_ID_2'], df_en_pro_rearranged['Gene_Name_2']\
    =np.vectorize(ReOrdering)(df_en_pro_small['Annotation_En_1'], df_en_pro_small['Gene_ID_1'], df_en_pro_small['Gene_Name_1'],\
    df_en_pro_small['Annotation_En_2'], df_en_pro_small['Gene_ID_2'], df_en_pro_small['Gene_Name_2'])

    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_en_pro_small.csv',df_en_pro_rearranged)
    #DECONVOLVE
    print('Analysing....')
    df1_c = df1.drop(['Gene_Name', 'TSS'], axis = 1)
    df2_c = df2.drop(['Gene_Name', 'TSS'], axis = 1)
    
    df1_decov1 = DeconvolveDataFrame(df1_c, 'Annotation_En', 'Code')
    print(df1_decov1.head())
    df1_decov2 = DeconvolveDataFrame(df1_decov1, 'Gene_ID', 'Code')
    
    df2_decov1 = DeconvolveDataFrame(df2_c, 'Annotation_En', 'Code')
    df2_decov2 = DeconvolveDataFrame(df2_decov1, 'Gene_ID', 'Code')

    df_decov_f = pd.merge(df1_decov2,df2_decov2, on = 'Code', how = 'outer', suffixes = ('_1','_2'))

    print(df_decov_f.head())

    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_Deconvolved.csv',df_decov_f)

    df_decov_small = df_decov_f[['File_Name_1','Code','Annotation_En_1','Gene_ID_1', 'Annotation_En_2','Gene_ID_2']]

    df_decov_small['Annotation_En_1_n'],df_decov_small['Gene_ID_1_n'],df_decov_small['Annotation_En_2_n'],df_decov_small['Gene_ID_2_n'] =\
    np.vectorize(ReorderingPromoter)(df_decov_small['Annotation_En_1'],\
        df_decov_small['Gene_ID_1'],df_decov_small['Annotation_En_2'],df_decov_small['Gene_ID_2'])

    df_c = df_decov_small.drop(['Annotation_En_1','Gene_ID_1', 'Annotation_En_2','Gene_ID_2'], axis = 1)

    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_Deconvolved_allign.csv',df_c)

    df_c_prom = df_c[df_c['Gene_ID_2_n'] != 'None']
    df_c_prom.drop_duplicates(inplace = True)

    df_c_prom_en = df_c_prom[df_c_prom['Annotation_En_1_n'] != 'None']
    df_c_prom_en = df_c_prom_en[['File_Name_1','Annotation_En_1_n','Gene_ID_2_n']]
    df_c_prom_en1 = df_c_prom_en.drop_duplicates()
    df_f = pd.DataFrame(df_c_prom_en1.groupby(['Annotation_En_1_n', 'Gene_ID_2_n']).count())
    df_f1 = df_f.reset_index(level =1)
    df_f2 = df_f1.reset_index()
    df_f2_f = df_f2[df_f2['File_Name_1'] > 1]
    #important_prom_en = df_f2_f.drop(['Gene_ID_1_n','Annotation_En_2_n'], axis = 1)
    
    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_Promoters_Enhancers_Significant.csv',df_f2_f)
    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_Deconvolved_Promoters.csv',df_c_prom)
    SaveDFtoCSV('Results_Super_Enhancers','Chia_Pet_Annotated_Deconvolved_Promoters_Enhancers.csv',df_c_prom_en)
    print('End')
    

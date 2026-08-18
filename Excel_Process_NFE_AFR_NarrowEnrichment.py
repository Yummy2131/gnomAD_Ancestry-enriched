import os
import pandas as pd

def filter(input_folder, output_file, group_column, enriched_frequency_column, B_frequency_column):
    # Initialize an empty list to collect DataFrames
    all_results = []

    #The left most column
    all_results.append(pd.DataFrame({'': ['FileName',
                                          '',
                                          'Total popmax #',
                                          'enriched #',
                                          '',
                                          'P/LP #',
                                          'B/LB #',
                                          'VUS #',
                                          'Conflicting #',
                                          'Blank #',
                                          'Other #',
                                          '',
                                          'pLOF',
                                          'CADD >/= 28.1',
                                          'CADD >/= 25.3',
                                          'CADD >/= 20',
                                          'REVEL >/= 0.932',
                                          'REVEL >/= 0.773',
                                          'REVEL >/= 0.644',
                                          'REVEL >/= 0.4',
                                          'SpliceAI >/= 0.8',
                                          'SpliceAI >/= 0.5',
                                          'SpliceAI >/= 0.2',
                                          '',
                                          'Enriched P/LP #',
                                          'Enriched B/LB #',
                                          'Enriched VUS #',
                                          'Enriched Conflicting #',
                                          'Enriched Blank #',
                                          'Enriched Other #',
                                          '',
                                          'Enriched pLOF',
                                          'Enriched CADD >/= 28.1',
                                          'Enriched CADD >/= 25.3',
                                          'Enriched CADD >/= 20',
                                          'Enriched REVEL >/= 0.932',
                                          'Enriched REVEL >/= 0.773',
                                          'Enriched REVEL >/= 0.644',
                                          'Enriched REVEL >/= 0.4',
                                          'Enriched SpliceAI >/= 0.8',
                                          'Enriched SpliceAI >/= 0.5',
                                          'Enriched SpliceAI >/= 0.2']}))

    # Loop through all files in the folder
    for filename in os.listdir(input_folder):
        if filename.endswith('.xls') or filename.endswith('.xlsx'):
            file_path = os.path.join(input_folder, filename)
            
            try:
                # Read the Excel file
                df = pd.read_excel(file_path)
                df.replace('#DIV/0!', pd.NA, inplace=True)
                #print(df)
                
                # GroupMax FAF group == 'Group_A' and 0.01 > Frequency > 0.00001
                filtered_A = df[(df[group_column] == Group_A) &
                                  (df[A_frequency_column] > 0.00005) &
                                  (df[A_frequency_column] < 0.01)]
                
                # GroupMax FAF group == 'Group_B' and 0.01 > Frequency > 0.00001
                filtered_B = df[(df[group_column] == Group_B) &
                                  (df[B_frequency_column] > 0.00005) &
                                  (df[B_frequency_column] < 0.01)]

                # GroupMax FAF group == 'Group_A' and 0.01 > Frequency > 0.00001 and Group_B AF / Group_A AF =≥10 or DIV/0! and Allele count ≥2
                Group_A_Enriched = df[
                                    (df[group_column] == Group_A) &
                                    (df[A_frequency_column] > 0.00005) &
                                    (df[A_frequency_column] < 0.01) &
                                    (df[A_allele_count] >= 2) &
                                    ((df[A_to_B_frequency_ratio_column] >= 10) | (df[A_to_B_frequency_ratio_column].isna())) &
                                    ((df[A_to_C_frequency_ratio_column] >= 10) | (df[A_to_C_frequency_ratio_column].isna()))]

                Group_B_Enriched = df[(df[group_column] == Group_B) &
                                      (df[B_frequency_column] > 0.00005) &
                                      (df[B_frequency_column] < 0.01) &
                                      (df[B_allele_count] >= 2) &
                                      ((df[B_to_A_frequency_ratio_column] >= 10) | (df[B_to_A_frequency_ratio_column].isna())) &
                                      ((df[B_to_C_frequency_ratio_column] >= 10) | (df[B_to_C_frequency_ratio_column].isna()))]

# -------------------------------------------------------------------------------------------------------------------------------------

                # Group_A popMax clinVar annotation
                Group_A_P_LP = filtered_A[filtered_A[clinVar_annotation].isin(P_LP)]

                Group_A_B_LB = filtered_A[filtered_A[clinVar_annotation].isin(B_LB)]

                Group_A_VUS = filtered_A[(filtered_A[clinVar_annotation] == 'Uncertain significance')]

                Group_A_Conflicting = filtered_A[(filtered_A[clinVar_annotation] == 'Conflicting classifications of pathogenicity')]

                Group_A_Blank = filtered_A[(filtered_A[clinVar_annotation].isna())]



                # Group_A popMax pLOF
                Group_A_pLOF = filtered_A[(filtered_A['VEP Annotation'] == 'frameshift_variant') |
                                          (filtered_A['VEP Annotation'] == 'stop_gained') |
                                          (filtered_A['VEP Annotation'] == 'startloss') |
                                          (filtered_A['VEP Annotation'] == 'stoploss') |
                                          (filtered_A['VEP Annotation'] == 'splice_acceptor_region') |
                                          (filtered_A['VEP Annotation'] == 'splice_donor_variant')]

                # Group_A CADD REVEL SpliceAI
                Group_A_CADD_28 = filtered_A[(filtered_A['cadd'] >= 28.1)]
                Group_A_CADD_25 = filtered_A[(filtered_A['cadd'] >= 25.3)]
                Group_A_CADD_20 = filtered_A[(filtered_A['cadd'] >= 20)]
                Group_A_REVEL_9 = filtered_A[(filtered_A['revel_max'] >= 0.932)]
                Group_A_REVEL_7 = filtered_A[(filtered_A['revel_max'] >= 0.773)]
                Group_A_REVEL_6 = filtered_A[(filtered_A['revel_max'] >= 0.644)]
                Group_A_REVEL_4 = filtered_A[(filtered_A['revel_max'] >= 0.4)]
                Group_A_SpliceAI_8 = filtered_A[(filtered_A['spliceai_ds_max'] >= 0.8)]
                Group_A_SpliceAI_5 = filtered_A[(filtered_A['spliceai_ds_max'] >= 0.5)]
                Group_A_SpliceAI_2 = filtered_A[(filtered_A['spliceai_ds_max'] >= 0.2)]

                # Group_B popMax ClinVar Annotation
                Group_B_P_LP = filtered_B[filtered_B[clinVar_annotation].isin(P_LP)]

                Group_B_B_LB = filtered_B[filtered_B[clinVar_annotation].isin(B_LB)]
                
                Group_B_VUS = filtered_B[(filtered_B[clinVar_annotation] == 'Uncertain significance')]

                Group_B_Conflicting = filtered_B[(filtered_B[clinVar_annotation] == 'Conflicting classifications of pathogenicity')]
                
                Group_B_Blank = filtered_B[(filtered_B[clinVar_annotation].isna())]


                # Group_B pLOF
                Group_B_pLOF = filtered_B[(filtered_B['VEP Annotation'] == 'frameshift_variant') |
                                        (filtered_B['VEP Annotation'] == 'stop_gained') |
                                        (filtered_B['VEP Annotation'] == 'startloss') |
                                        (filtered_B['VEP Annotation'] == 'stoploss') |
                                        (filtered_B['VEP Annotation'] == 'splice_acceptor_region') |
                                        (filtered_B['VEP Annotation'] == 'splice_donor_variant')]

                # Group_B popMax CADD REVEL SpliceAI
                Group_B_CADD_28 = filtered_B[(filtered_B['cadd'] >= 28.1)]
                Group_B_CADD_25 = filtered_B[(filtered_B['cadd'] >= 25.3)]
                Group_B_CADD_20 = filtered_B[(filtered_B['cadd'] >= 20)]
                Group_B_REVEL_9 = filtered_B[(filtered_B['revel_max'] >= 0.932)]
                Group_B_REVEL_7 = filtered_B[(filtered_B['revel_max'] >= 0.773)]
                Group_B_REVEL_6 = filtered_B[(filtered_B['revel_max'] >= 0.644)]
                Group_B_REVEL_4 = filtered_B[(filtered_B['revel_max'] >= 0.4)]
                Group_B_SpliceAI_8 = filtered_B[(filtered_B['spliceai_ds_max'] >= 0.8)]
                Group_B_SpliceAI_5 = filtered_B[(filtered_B['spliceai_ds_max'] >= 0.5)]
                Group_B_SpliceAI_2 = filtered_B[(filtered_B['spliceai_ds_max'] >= 0.2)]

#-------------------------------------------------------------------------------------------------------------------------------------
                # Enriched Group A ClinVar Annotation
                Group_A_P_LP_Enriched = Group_A_Enriched[Group_A_Enriched[clinVar_annotation].isin(P_LP)]

                Group_A_B_LB_Enriched = Group_A_Enriched[Group_A_Enriched[clinVar_annotation].isin(B_LB)]

                Group_A_VUS_Enriched = Group_A_Enriched[(Group_A_Enriched[clinVar_annotation] == 'Uncertain significance')]

                Group_A_Conflicting_Enriched = Group_A_Enriched[
                    (Group_A_Enriched[clinVar_annotation] == 'Conflicting classifications of pathogenicity')]

                Group_A_Blank_Enriched = Group_A_Enriched[(Group_A_Enriched[clinVar_annotation].isna())]


                # Enriched Group_A pLOF
                Group_A_pLOF_Enriched = Group_A_Enriched[(Group_A_Enriched['VEP Annotation'] == 'frameshift_variant') |
                                                         (Group_A_Enriched['VEP Annotation'] == 'stop_gained') |
                                                         (Group_A_Enriched['VEP Annotation'] == 'startloss') |
                                                         (Group_A_Enriched['VEP Annotation'] == 'stoploss') |
                                                         (Group_A_Enriched['VEP Annotation'] == 'splice_acceptor_region') |
                                                         (Group_A_Enriched['VEP Annotation'] == 'splice_donor_variant')]

                # Enriched Group_A CADD REVEL SpliceAI
                Group_A_CADD_28_Enriched = Group_A_Enriched[(Group_A_Enriched['cadd'] >= 28.1)]
                Group_A_CADD_25_Enriched = Group_A_Enriched[(Group_A_Enriched['cadd'] >= 25.3)]
                Group_A_CADD_20_Enriched = Group_A_Enriched[(Group_A_Enriched['cadd'] >= 20)]
                Group_A_REVEL_9_Enriched = Group_A_Enriched[(Group_A_Enriched['revel_max'] >= 0.932)]
                Group_A_REVEL_7_Enriched = Group_A_Enriched[(Group_A_Enriched['revel_max'] >= 0.773)]
                Group_A_REVEL_6_Enriched = Group_A_Enriched[(Group_A_Enriched['revel_max'] >= 0.644)]
                Group_A_REVEL_4_Enriched = Group_A_Enriched[(Group_A_Enriched['revel_max'] >= 0.4)]
                Group_A_SpliceAI_8_Enriched = Group_A_Enriched[(Group_A_Enriched['spliceai_ds_max'] >= 0.8)]
                Group_A_SpliceAI_5_Enriched = Group_A_Enriched[(Group_A_Enriched['spliceai_ds_max'] >= 0.5)]
                Group_A_SpliceAI_2_Enriched = Group_A_Enriched[(Group_A_Enriched['spliceai_ds_max'] >= 0.2)]

                # Enriched Group B ClinVar Annotation
                Group_B_P_LP_Enriched = Group_B_Enriched[Group_B_Enriched[clinVar_annotation].isin(P_LP)]

                Group_B_B_LB_Enriched = Group_B_Enriched[Group_B_Enriched[clinVar_annotation].isin(B_LB)]

                Group_B_VUS_Enriched = Group_B_Enriched[(Group_B_Enriched[clinVar_annotation] == 'Uncertain significance')]

                Group_B_Conflicting_Enriched = Group_B_Enriched[(Group_B_Enriched[clinVar_annotation] == 'Conflicting classifications of pathogenicity')]

                Group_B_Blank_Enriched = Group_B_Enriched[(Group_B_Enriched[clinVar_annotation].isna())]


                # Enriched Group_B pLOF
                Group_B_pLOF_Enriched = Group_B_Enriched[(Group_B_Enriched['VEP Annotation'] == 'frameshift_variant') |
                                                         (Group_B_Enriched['VEP Annotation'] == 'stop_gained') |
                                                         (Group_B_Enriched['VEP Annotation'] == 'startloss') |
                                                         (Group_B_Enriched['VEP Annotation'] == 'stoploss') |
                                                         (Group_B_Enriched['VEP Annotation'] == 'splice_acceptor_region') |
                                                         (Group_B_Enriched['VEP Annotation'] == 'splice_donor_variant')]

                # Enriched Group_B CADD REVEL SpliceAI
                Group_B_CADD_28_Enriched = Group_B_Enriched[(Group_B_Enriched['cadd'] >= 28.1)]
                Group_B_CADD_25_Enriched = Group_B_Enriched[(Group_B_Enriched['cadd'] >= 25.3)]
                Group_B_CADD_20_Enriched = Group_B_Enriched[(Group_B_Enriched['cadd'] >= 20)]
                Group_B_REVEL_9_Enriched = Group_B_Enriched[(Group_B_Enriched['revel_max'] >= 0.932)]
                Group_B_REVEL_7_Enriched = Group_B_Enriched[(Group_B_Enriched['revel_max'] >= 0.773)]
                Group_B_REVEL_6_Enriched = Group_B_Enriched[(Group_B_Enriched['revel_max'] >= 0.644)]
                Group_B_REVEL_4_Enriched = Group_B_Enriched[(Group_B_Enriched['revel_max'] >= 0.4)]
                Group_B_SpliceAI_8_Enriched = Group_B_Enriched[(Group_B_Enriched['spliceai_ds_max'] >= 0.8)]
                Group_B_SpliceAI_5_Enriched = Group_B_Enriched[(Group_B_Enriched['spliceai_ds_max'] >= 0.5)]
                Group_B_SpliceAI_2_Enriched = Group_B_Enriched[(Group_B_Enriched['spliceai_ds_max'] >= 0.2)]
                
                # Count the number of rows that match each condition
                Group_A_row_count = len(filtered_A)
                Group_B_row_count = len(filtered_B)
                Group_A_Enriched_count = len(Group_A_Enriched)
                Group_B_Enriched_count = len(Group_B_Enriched)

                # Group_A popMax Count
                Group_A_pLOF_count = len(Group_A_pLOF)
                Group_A_P_LP_count = len(Group_A_P_LP)
                Group_A_B_LB_count = len(Group_A_B_LB)
                Group_A_VUS_count = len(Group_A_VUS)
                Group_A_Conflicting_count = len(Group_A_Conflicting)
                Group_A_Blank_count = len(Group_A_Blank)
                Group_A_Other_count = Group_A_row_count-Group_A_P_LP_count-Group_A_B_LB_count-Group_A_VUS_count-Group_A_Conflicting_count-Group_A_Blank_count
                Group_A_CADD_28_count = len(Group_A_CADD_28)
                Group_A_CADD_25_count = len(Group_A_CADD_25)
                Group_A_CADD_20_count = len(Group_A_CADD_20)
                Group_A_REVEL_9_count = len(Group_A_REVEL_9)
                Group_A_REVEL_7_count = len(Group_A_REVEL_7)
                Group_A_REVEL_6_count = len(Group_A_REVEL_6)
                Group_A_REVEL_4_count = len(Group_A_REVEL_4)
                Group_A_SpliceAI_8_count = len(Group_A_SpliceAI_8)
                Group_A_SpliceAI_5_count = len(Group_A_SpliceAI_5)
                Group_A_SpliceAI_2_count = len(Group_A_SpliceAI_2)

                # Group_A Enriched Count
                Group_A_Enriched_pLOF_count = len(Group_A_pLOF_Enriched)
                Group_A_Enriched_P_LP_count = len(Group_A_P_LP_Enriched)
                Group_A_Enriched_B_LB_count = len(Group_A_B_LB_Enriched)
                Group_A_Enriched_VUS_count = len(Group_A_VUS_Enriched)
                Group_A_Enriched_Conflicting_count = len(Group_A_Conflicting_Enriched)
                Group_A_Enriched_Blank_count = len(Group_A_Blank_Enriched)
                Group_A_Enriched_Other_count = Group_A_Enriched_count-Group_A_Enriched_P_LP_count-Group_A_Enriched_B_LB_count-Group_A_Enriched_VUS_count-Group_A_Enriched_Conflicting_count-Group_A_Enriched_Blank_count
                Group_A_Enriched_CADD_28_count = len(Group_A_CADD_28_Enriched)
                Group_A_Enriched_CADD_25_count = len(Group_A_CADD_25_Enriched)
                Group_A_Enriched_CADD_20_count = len(Group_A_CADD_20_Enriched)
                Group_A_Enriched_REVEL_9_count = len(Group_A_REVEL_9_Enriched)
                Group_A_Enriched_REVEL_7_count = len(Group_A_REVEL_7_Enriched)
                Group_A_Enriched_REVEL_6_count = len(Group_A_REVEL_6_Enriched)
                Group_A_Enriched_REVEL_4_count = len(Group_A_REVEL_4_Enriched)
                Group_A_Enriched_SpliceAI_8_count = len(Group_A_SpliceAI_8_Enriched)
                Group_A_Enriched_SpliceAI_5_count = len(Group_A_SpliceAI_5_Enriched)
                Group_A_Enriched_SpliceAI_2_count = len(Group_A_SpliceAI_2_Enriched)

                
                # Group_B popMax Count
                Group_B_pLOF_count = len(Group_B_pLOF)
                Group_B_P_LP_count = len(Group_B_P_LP)
                Group_B_B_LB_count = len(Group_B_B_LB)
                Group_B_VUS_count = len(Group_B_VUS)
                Group_B_Conflicting_count = len(Group_B_Conflicting)
                Group_B_Blank_count = len(Group_B_Blank)
                Group_B_Other_count = Group_B_row_count-Group_B_P_LP_count-Group_B_B_LB_count-Group_B_VUS_count-Group_B_Conflicting_count-Group_B_Blank_count
                Group_B_CADD_28_count = len(Group_B_CADD_28)
                Group_B_CADD_25_count = len(Group_B_CADD_25)
                Group_B_CADD_20_count = len(Group_B_CADD_20)
                Group_B_REVEL_9_count = len(Group_B_REVEL_9)
                Group_B_REVEL_7_count = len(Group_B_REVEL_7)
                Group_B_REVEL_6_count = len(Group_B_REVEL_6)
                Group_B_REVEL_4_count = len(Group_B_REVEL_4)
                Group_B_SpliceAI_8_count = len(Group_B_SpliceAI_8)
                Group_B_SpliceAI_5_count = len(Group_B_SpliceAI_5)
                Group_B_SpliceAI_2_count = len(Group_B_SpliceAI_2)

                # Group_B Enriched Count
                Group_B_Enriched_pLOF_count = len(Group_B_pLOF_Enriched)
                Group_B_Enriched_P_LP_count = len(Group_B_P_LP_Enriched)
                Group_B_Enriched_B_LB_count = len(Group_B_B_LB_Enriched)
                Group_B_Enriched_VUS_count = len(Group_B_VUS_Enriched)
                Group_B_Enriched_Conflicting_count = len(Group_B_Conflicting_Enriched)
                Group_B_Enriched_Blank_count = len(Group_B_Blank_Enriched)
                Group_B_Enriched_Other_count = Group_B_Enriched_count-Group_B_Enriched_P_LP_count-Group_B_Enriched_B_LB_count-Group_B_Enriched_VUS_count-Group_B_Enriched_Conflicting_count-Group_B_Enriched_Blank_count
                Group_B_Enriched_CADD_28_count = len(Group_B_CADD_28_Enriched)
                Group_B_Enriched_CADD_25_count = len(Group_B_CADD_25_Enriched)
                Group_B_Enriched_CADD_20_count = len(Group_B_CADD_20_Enriched)
                Group_B_Enriched_REVEL_9_count = len(Group_B_REVEL_9_Enriched)
                Group_B_Enriched_REVEL_7_count = len(Group_B_REVEL_7_Enriched)
                Group_B_Enriched_REVEL_6_count = len(Group_B_REVEL_6_Enriched)
                Group_B_Enriched_REVEL_4_count = len(Group_B_REVEL_4_Enriched)
                Group_B_Enriched_SpliceAI_8_count = len(Group_B_SpliceAI_8_Enriched)
                Group_B_Enriched_SpliceAI_5_count = len(Group_B_SpliceAI_5_Enriched)
                Group_B_Enriched_SpliceAI_2_count = len(Group_B_SpliceAI_2_Enriched)
                
                # Create a DataFrame with the file name, headers, counts
                result_df = pd.DataFrame({
                    f'{filename}': [filename,
                                    Group_A,
                                    Group_A_row_count,
                                    Group_A_Enriched_count,
                                    '',
                                    Group_A_P_LP_count,
                                    Group_A_B_LB_count,
                                    Group_A_VUS_count,
                                    Group_A_Conflicting_count,
                                    Group_A_Blank_count,
                                    Group_A_Other_count,
                                    '',
                                    Group_A_pLOF_count,
                                    Group_A_CADD_28_count,
                                    Group_A_CADD_25_count,
                                    Group_A_CADD_20_count,
                                    Group_A_REVEL_9_count,
                                    Group_A_REVEL_7_count,
                                    Group_A_REVEL_6_count,
                                    Group_A_REVEL_4_count,
                                    Group_A_SpliceAI_8_count,
                                    Group_A_SpliceAI_5_count,
                                    Group_A_SpliceAI_2_count,
                                    '',
                                    Group_A_Enriched_P_LP_count,
                                    Group_A_Enriched_B_LB_count,
                                    Group_A_Enriched_VUS_count,
                                    Group_A_Enriched_Conflicting_count,
                                    Group_A_Enriched_Blank_count,
                                    Group_A_Enriched_Other_count,
                                    '',
                                    Group_A_Enriched_pLOF_count,
                                    Group_A_Enriched_CADD_28_count,
                                    Group_A_Enriched_CADD_25_count,
                                    Group_A_Enriched_CADD_20_count,
                                    Group_A_Enriched_REVEL_9_count,
                                    Group_A_Enriched_REVEL_7_count,
                                    Group_A_Enriched_REVEL_6_count,
                                    Group_A_Enriched_REVEL_4_count,
                                    Group_A_Enriched_SpliceAI_8_count,
                                    Group_A_Enriched_SpliceAI_5_count,
                                    Group_A_Enriched_SpliceAI_2_count],
                    
                    '': ['',
                         Group_B,
                         Group_B_row_count,
                         Group_B_Enriched_count,
                         '',
                         Group_B_P_LP_count,
                         Group_B_B_LB_count,
                         Group_B_VUS_count,
                         Group_B_Conflicting_count,
                         Group_B_Blank_count,
                         Group_B_Other_count,
                         '',
                         Group_B_pLOF_count,
                         Group_B_CADD_28_count,
                         Group_B_CADD_25_count,
                         Group_B_CADD_20_count,
                         Group_B_REVEL_9_count,
                         Group_B_REVEL_7_count,
                         Group_B_REVEL_6_count,
                         Group_B_REVEL_4_count,
                         Group_B_SpliceAI_8_count,
                         Group_B_SpliceAI_5_count,
                         Group_B_SpliceAI_2_count,
                         '',
                         Group_B_Enriched_P_LP_count,
                         Group_B_Enriched_B_LB_count,
                         Group_B_Enriched_VUS_count,
                         Group_B_Enriched_Conflicting_count,
                         Group_B_Enriched_Blank_count,
                         Group_B_Enriched_Other_count,
                         '',
                         Group_B_Enriched_pLOF_count,
                         Group_B_Enriched_CADD_28_count,
                         Group_B_Enriched_CADD_25_count,
                         Group_B_Enriched_CADD_20_count,
                         Group_B_Enriched_REVEL_9_count,
                         Group_B_Enriched_REVEL_7_count,
                         Group_B_Enriched_REVEL_6_count,
                         Group_B_Enriched_REVEL_4_count,
                         Group_B_Enriched_SpliceAI_8_count,
                         Group_B_Enriched_SpliceAI_5_count,
                         Group_B_Enriched_SpliceAI_2_count]
                })

                # Append the result dataframe to the list
                all_results.append(result_df)

                # Add an empty column between files by appending an empty DataFrame
                all_results.append(pd.DataFrame({'': ['', '', '', '', '', '', '', '', '', '', '', '']})) 

                print(f"{filename}: {Group_A_row_count} rows {Group_A}, {Group_B_row_count} rows {Group_B} match the conditions")
                
            except Exception as e:
                print(f"Error processing {filename}: {e}")
    
    # Concatenate all the result dataframes horizontally (axis=1) with empty columns between
    final_result_df = pd.concat(all_results, axis=1)

    # Write the result DataFrame to a new Excel file
    final_result_df.to_excel(output_file, index=False, header=False)
    print(f"Filtered saved to {output_file}")

# Variables
file_path = 'All Aggregated'
input_folder = f'C:\\Users\\boxx_\\Desktop\\gnomAD\\{file_path}\\xlsx'    # Folder containing the input Excel files
output_folder = f'C:\\Users\\boxx_\\Desktop\\gnomAD\\{file_path}'  # Folder where the filtered Excel files will be saved
output_file = f'C:\\Users\\boxx_\\Desktop\\gnomAD\\{file_path}\\NFE_AFR_{file_path}_NarrowEnrichment.xlsx'
group_column = 'GroupMax FAF group'
Group_A = 'nfe' #nfe or amr or afr
Group_B = 'afr' #nfe or amr or afr
A_frequency_column = 'Frequency European non-Finnish'  # Column for the A group frequency
B_frequency_column = 'Frequency African/African American'    #
A_allele_count = 'Allele Count European (non-Finnish)'
B_allele_count = 'Allele Count African/African American' #
A_to_B_frequency_ratio_column = 'European non-Finnish Frequency/African African American Frequency'
A_to_C_frequency_ratio_column = 'European non-Finnish Frequency/Admixed Frequency'
B_to_A_frequency_ratio_column =  'African African American Frequency/European non-Finnish Frequency'
B_to_C_frequency_ratio_column =  'African African American Frequency/Admixed Frequency'
clinVar_annotation = 'ClinVar Germline Classification' #ClinVar Germline Classification /ClinVar Clinical Significance
P_LP = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']
B_LB = ['Benign','Likely benign','Benign/Likely benign']
VUS = ['Uncertain significance']
CON =['Conflicting classifications of pathogenicity']


filter(input_folder, output_file, group_column, A_frequency_column, B_frequency_column)


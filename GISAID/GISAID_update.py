#!/usr/bin/env python

'''by James C. Hu
Change log

01.01.2023
    "New-year-new-me" update:
    - Merged ASU dashboard Exe scripts.
        - Functionality of each script converted to custom functions.

01.04.2023
    General:
        Script updated with new VOC and VBM list/dictionaries. Now reflect CDC covid tracker.
        (https://covid.cdc.gov/covid-data-tracker/#variant-proportions)

    Compare with Exe_GISAID_Update_v0.2.5.py to see full change.

    1)  Lineage lists and dictionaries:
            - Added new list:(U001)
                - Recombinant_lineages = ['XBB']
            - Updated list/dict: (U002)
                - Lists
                    from:
                    - All_lineages = ... Numeric_lineages
                    - All_lineages_list = [... Numeric_lineages]
                    to:
                    - All_lineages = ... Numeric_lineages + Recombinant_lineages
                    - All_lineages_list = [... Numeric_lineages, Recombinant_lineages]
                - Dictionaries
                    from:
                    - VOC_dict = {... 'BA.5':'Omicron (BA.5)'}
                    to:
                    - VOC_dict = {... 'BA.5':'Omicron (BA.5)', BQ.1.1': 'Omicron_BA.5.3.1.1.1.1.1.1 (BQ.1.1)', 'BQ.1': 'Omicron_BA.5.3.1.1.1.1.1 (BQ.1)',
                    'XBB': 'Recombinant_BA.2 (XBB)', 'BN.1': 'Omicron_BA.2.75.51 (BN.1)', 'BF.7': 'Omicron_BA.5.2.11 (BF.7)',
                    'BA.5.2.6': 'Omicron (BA.5.2.6)', 'BA.4.6': 'Omicron (BA.4.6)', 'BF.11': 'Omicron_BA.5.2.1.11 (BF.11)',
                    'BA.2.75.2': 'Omicron (BA.2.75.2)', 'BA.1.1': 'Omicron (BA.1.1)', 'B.1.1.529': 'Omicron (Omicron)'}

    2)  Function exicution:
            - Created a new df to house lineage_count function output for the newly created Recombinant_lineages list
                in the Lineage lists and dictionaries section. (U003)
                - df18 = pd.DataFrame(lineage_count(Recombinant_lineages, sorted_t2a_data_df))
            - Update df_final with the newly added sequential df.
                from:
                - df_final = pd.concat([df02, df03, df04, df05, df06, df07, df08,
                      df09, df10, df11, df12, df13, df14, df15, df16, df17], axis=1)
                to:
                - df_final = pd.concat([df02, df03, df04, df05, df06, df07, df08,
                      df09, df10, df11, df12, df13, df14, df15, df16, df17, df18], axis=1)

01.05.2023
    New mutation counting function added to Counting Function section @ HHHH
    Added code to write output log file.
    Added new Utitlity function to change color of terminal output.

01.06.2023
    New mutation counting funciton added to Counting Function section @ IIII

02.21.2023
    New VOC list added. Option to change VOC output @ JJJJ

08.02.2023
    due to pandas update, pd.to_datetime method now requires format argument.
    added format='mixed' to gisaid_data_sort_t1 function at line 135, gisaid_date_summary_count funciton at line 195

General Notes:
    IDE: Sublime 4
        - Package Control
        - BracketHighligher
        - Material Theme
        - Side Bar
        - Anaconda (Same name, but not the same functionality as the package control/distributor)
            - Features include auto pep8 formatting, function autofill, error prediciton, syntax suggestion

    Build System:  Python 3.9.15+
                    - Pandas 1.5.2+
                    - numpy 1.22.3+

Table of contents

    AAAA: File merging function
        - AAA1: Renaming merging, input, and output files.
    BBBB: Tier 1 parsing function
        - BBB1: List(s) for tier 1 parsing.
    CCCC: Tier 2 parsing funciton
        - CCC1: Lists for tier 2 parsing: Counting funciton for lineages by week.
    DDDD: Counting function for dates. Currently written to output earliest and most recent dates for input criteria.
    EEEE: Counting funciton for input. Currenlty written to count Lineages.
    FFFF: Utility function for dropping columns from a df.
    - FFF1: List as index for columns to be dropped.
    GGGG: Utility function for parsing rows from a df.
    - Currently being used to create summary tables by slicing from existing df.
    - GGG1: Lists and dictionaries used for summary tables.
    HHHH: The function below parses a tsv file for specified values within a given list
    IIII: Utility function for changing terminal output text color
    - III1: Change value to change output color.
    JJJJ: Select which VOC list to output.
'''
#######################################################################################################################
# Dependencies
#######################################################################################################################
import pandas as pd
import numpy as np
import sys
from datetime import datetime
from collections import Counter
pd.options.display.width = 0  # Displays full terminal output
pd.options.mode.chained_assignment = None
# Stops the chain_assignemnt error from popping up.
# Currently does not contribute to misscount, but I should look into it more deeply at some point.

#######################################################################################################################
# Parsing Functions
#######################################################################################################################

'''
The function below performs tier 1 parsing on selected tsv file.
Tier 1 parsing includes:
    - Slicing out desired columns.
    - Converting date formats.
    note: all hard-coded variables can be made more dynamic.
outfile1 = new combined storage file
t1_columns = desired columns for t1 parsing
BBBB
'''


def gisaid_data_sort_t1(outfile1, t1_columns):
    linage_list = []
    df = pd.read_csv(outfile1, sep='\t', header=0, dtype=str)
    linage_list = [i for i in df['Lineage'] if i not in linage_list]
    df1 = df[t1_columns].copy().set_index('Virus name').fillna('NA')
    df1['Lineage'] = df1['Lineage'].str.replace('None', 'NA')
    df1['Date'] = pd.to_datetime(df1['Collection date'], format='mixed')
    df1['Ouptput_date_format'] = [datetime.strptime(
        str(x), '%Y-%m-%d').strftime('%b-%d-%Y') for x in df1['Date'].dt.date]
    df1['Week_number'] = df1['Date'].dt.isocalendar().week
    df1['Year'] = df1['Date'].dt.year
    df1['First_day_of_Week'] = df1['Date'] - \
        df1['Date'].dt.weekday * np.timedelta64(1, 'D')
    df1['Week'] = [datetime.strptime(str(x), '%Y-%m-%d').strftime('%b-%d-%Y')
                   for x in df1['First_day_of_Week'].dt.date]
    df1.sort_values(by='First_day_of_Week', ascending=True, inplace=True)
    return df1


'''
The function below performs tier 2 parsing on selected tsv file.
Tier 2 parsing includes:
    - generating counts for selected column by time value.
    - Converting date formats.
    note: all hard-coaded variables can be made more dynamic.
t2a_columns = desired columns for t2a parsing
count_column_str = column name used for counting, as a string
Date_or_Week_str = Count by 'Date' or 'Week' depending on which column exists in the df.
CCCC
'''


def gisaid_data_sort_t2(sorted_t1_data_df, t2a_columns, count_column_str, Date_or_Week_str):
    df = sorted_t1_data_df.groupby(t2a_columns)[
        count_column_str].count().reset_index(name='Count')
    df[Date_or_Week_str] = pd.to_datetime(df[Date_or_Week_str])
    df = df.sort_values(by=Date_or_Week_str).reset_index(drop=True).pivot(
        index=Date_or_Week_str, columns=count_column_str, values='Count').fillna(0)
    df.columns = [i for i in df]
    df.reset_index(inplace=True)
    df[Date_or_Week_str] = pd.to_datetime(df[Date_or_Week_str])
    df[Date_or_Week_str] = [datetime.strptime(
        str(x), '%Y-%m-%d').strftime('%b-%d-%Y') for x in df[Date_or_Week_str].dt.date]
    df.set_index(Date_or_Week_str, inplace=True)
    return df


#######################################################################################################################
# Counting Functions
#######################################################################################################################

'''
The function below cleans the tier 2 output; includes count.
    - Displays the indexing column for counts, the count, as well as formated earlist and most recent dates.
note: funciton input is currently named sorted_t2b_data_df for clarity; input can be any df, but code will need
to be adjusted accordingly
DDDD
'''


def gisaid_date_summary_count(sorted_t2b_data_df):
    df = pd.DataFrame(sorted_t2b_data_df.sum().transpose()).reset_index()
    df.rename(columns={'index': 'Lineage',
                       0: 'Cumulative_count'}, inplace=True)
    df['Earliest_case'], df['Most_recent_case'] = 'NA', 'NA'
    for i in sorted_t2b_data_df:
        dfa = pd.to_datetime(
            sorted_t2b_data_df[sorted_t2b_data_df[i] > 0].index.tolist(), format='mixed')
        date_list.append([i, dfa.min(), dfa.max()])
    for i in date_list:
        df['Earliest_case'].loc[(df['Lineage'] == i[0])] = i[1]
        df['Most_recent_case'].loc[(df['Lineage'] == i[0])] = i[2]
    df['Earliest_case'] = pd.to_datetime(df['Earliest_case'])
    df['Most_recent_case'] = pd.to_datetime(df['Most_recent_case'])
    df['Earliest_case'] = [datetime.strptime(
        str(x), '%Y-%m-%d').strftime('%b-%d-%Y') for x in df['Earliest_case'].dt.date]
    df['Most_recent_case'] = [datetime.strptime(
        str(x), '%Y-%m-%d').strftime('%b-%d-%Y') for x in df['Most_recent_case'].dt.date]
    date_list.clear()
    df.set_index('Lineage', inplace=True)
    return df


'''
The function below performs the counting for each column of the sorted tier 2 data.
specific_v_list = Alpha_lineages, Beta_lineages, Gamma_lineages etc
EEEE
Counting funciton for counting input criteria. Currenlty written to count Lneages.
'''


def lineage_count(specific_lineage_list, sorted_t2a_data_df):
    dfs = []
    df = pd.DataFrame()
    for target_lineage in specific_lineage_list:
        try:
            if specific_lineage_list[0] not in sorted_t2a_data_df and df[specific_lineage_list[0]] not in dfs:
                df = pd.DataFrame(
                    columns=[specific_lineage_list[0]], index=sorted_t2a_data_df.index).fillna(0)
                dfs.append(df)
            else:
                dfs.append(sorted_t2a_data_df[[specific_lineage_list[0]]])
        except Exception:
            pass
    for target_lineage in specific_lineage_list:
        try:
            if specific_lineage_list[0] not in sorted_t2a_data_df:
                df = pd.DataFrame(
                    columns=[specific_lineage_list[0]], index=sorted_t2a_data_df.index).fillna(0)
                dfs.append(df)
            else:
                dfs.append(sorted_t2a_data_df[[specific_lineage_list[0]]])
        except Exception:
            pass
    for target_lineage in specific_lineage_list:
        # print(target_lineage) # print should return the lineages within the specified list
        for lineage in sorted_t2a_data_df:
            # print(lineage) # print should return the names of the columns availble for parsing.
            try:
                if str(lineage).startswith(target_lineage):
                    # print(lineage, lineage_check) # Unhash to see match behavior
                    dfs.append(sorted_t2a_data_df[[lineage]])
            except Exception:
                pass
    dfs = pd.concat(dfs, axis=1)
    # dfs.to_csv('dfs_before.csv')
    dfs = dfs.loc[:, ~dfs.columns.duplicated()]
    # dfs.to_csv('dfs_after.csv')
    sum_column = dfs.sum(axis=1).reset_index(
        name=specific_lineage_list[0]).set_index('Week')
    return sum_column


'''
HHHH
The function below parses a tsv file for specified values within a given list.
Currently coded to parse updated combined file AA Substitutions column for mutations.
'''


def mutation_counter(tsv, primary_mutation_list):
    filtered_mutation_list = []
    df = pd.read_csv(tsv, sep='\t', header=0, dtype=str)
    total_sample_count = len(df)
    for mutation_csv in df['AA Substitutions']:
        mutation_csv = mutation_csv.strip('()')
        mutation_csv = mutation_csv.split(',')
        # print(mutation_csv)
        for mutation in mutation_csv:
            for target_mutation in primary_mutation_list:
                if mutation.startswith(target_mutation):
                    filtered_mutation_list.append(mutation)
    mutation_and_counts = Counter(filtered_mutation_list)
    filtered_mutation_set = set(filtered_mutation_list)
    filtered_mutation_list = list(filtered_mutation_set)
    filtered_mutation_list.sort()
    total_mutation_count = sum(mutation_and_counts.values())
    df1 = pd.DataFrame(columns=['Total_sample_count', 'Total_muatations_at_codon', 'Specific_muation_count',
                                'Specific_mutation_prevalence', 'Total_prevalence'], index=filtered_mutation_list)
    for i, j in zip(filtered_mutation_list, mutation_and_counts.values()):
        df1['Total_sample_count'].iloc[(df1.index == i)] = total_sample_count
        df1['Total_muatations_at_codon'].iloc[(
            df1.index == i)] = total_mutation_count
        df1['Specific_muation_count'].iloc[(df1.index == i)] = j
        df1['Specific_mutation_prevalence'].iloc[(
            df1.index == i)] = f'{round(((j / total_mutation_count) * 100), 2)} %'
        df1['Total_prevalence'].iloc[(
            df1.index == i)] = f'{round(((j / total_sample_count) * 100), 2)} %'
    df1.to_csv(f'{outdate}_Mutation_summary.csv')
    print('\nPrinting: mutation_counter function output')
    return print(df1)


#######################################################################################################################
# Utility Functions
#######################################################################################################################


'''
AAAA
The function below merges two csv files by index.
Index is currently hard-coaded to be 'Virus name' but can be changed to be more dynamic if needed.
infile1 = combined storage file
infile2 = new data file
'''


def gisaid_update(infile1, infile2):
    df = pd.read_csv(infile1, sep='\t', header=0, dtype=str)
    df = df.set_index('Virus name')
    df['Additional location information'] = df['Additional location information'].astype(
        str)
    df['Unnamed: 13'] = df['Unnamed: 13'].astype(str)
    df1 = pd.read_csv(infile2, sep='\t', header=0, dtype=str)
    df1 = df1.set_index('Virus name')
    # df1['Additional location information'] = df1['Additional location information'].astype(
    #     str)
    df['Unnamed: 13'] = df['Unnamed: 13'].astype(str)
    df2 = pd.concat([df, df1])
    df2 = df2.drop_duplicates()
    return df2


'''
FFFF
The function below dropkicks undesierable df columns out of existance.
All_lineages_list = All_lineages
'''


def column_drop(df, All_lineages_list):
    for i in All_lineages_list:
        for j in df:
            if j.startswith(i):
                try:
                    df = df.drop(columns=[j])
                except Exception:
                    pass
    return df


'''
GGGG
The function below parses column values by index from selected dataframe.
variant_of_list = list for VOC, VOI, or VBM
variant_of_dict = dict for VOC, VOI, or VBM
'''


def table_summary(df, variant_of_list, variant_of_dict):
    df1 = pd.DataFrame(
        columns=['Lineage', 'Cumulative_count', 'Earliest_case', 'Most_recent_case']).set_index('Lineage')
    for i in variant_of_list:
        try:
            df1.loc[i] = df.loc[i]
        except Exception:
            pass
    df1.rename(index=dict(variant_of_dict), inplace=True)
    return df1


'''
IIII
The function below changes the color of a string output to the terminal window.
'''


def color_me(str_variable, number):
    # color options: 7, 30-39, 91, 92
    color = f'\033[{number}m'
    end = '\033[0m'
    return print(color + str_variable + end)


#######################################################################################################################
# Variables and funciton inputs
#######################################################################################################################


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
# File Naming
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
outdate = datetime.today().strftime('%Y%m%d')

stdoutOrigin = sys.stdout
sys.stdout = open(f'{outdate}_GISAID_update_log.txt', 'w')

# AAA1
# Current GISAID records
infile1 = 'CombinedGISAID_hcov-19_2023_08_22_16.tsv'
# Latest GISAID records
infile2 = 'gisaid_hcov-19_2023_08_29_16.tsv'
# Updated Gisaid records
outfile1 = 'CombinedGISAID_hcov-19_2023_08_29_16.tsv'
outfile2 = f'{outdate}_lineage_by_week_report.csv'
outfile3 = f'{outdate}_GISAID_date_summary.csv'
# Unhash for GISAID t1 intermediate file
outfile4 = f'{outdate}_GISAID_t1_output.csv'
# outfile5 = f'{outdate}_GISAID_t2_output.csv' # Unhash for GISAID t2 intermediate file


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
# Sorting function inputs
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
'''
Notes:
    - New queries can exist as new lists housing desired columns to be parsed.
    e.g. tier 2 sorting of 'Lineage' count by 'Week' or 'Output_date_format'
'''
t1_columns = ['Virus name', 'Accession ID',
              'Collection date', 'Lineage', 'AA Substitutions']  # BBB1

t2a_columns = ['Lineage', 'Week']  # CCC1
t2b_columns = ['Lineage', 'Ouptput_date_format']

linage_list = []  # used for creating lineage_by_week df.
drop_list = []  # used to drop specific linages from count.
date_list = []  # used for combining datetime objects for weekly sort.


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
# Lineage lists and dictionaries
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
'''
EEEE
Add new PANGO designations as new lists below, then update All_lineages and All_lineages_list to include new list.
i.e. NewPangoDesignation_lineages = []
Add new sublineages to their respective parent lineage.
i.e. NewPangoDesignation_lineages = ['NewPangoLineage', 'NewPangoSubLineage']

!!!NOTE!!!  For best results, list[0] should to be the parent (B.1.1.7, B.1.351, etc.) lineage !!!NOTE!!!

Add variables into their respecive lists in ways that are compatible with .startswith().
e.g. Adding the . in Q. for alpha lineages.

Error Example: B.1.617 and B.1.617.1 are != and both exist somewhere within a lineage list.
               B.1.617 will include the count for B.1.617.1
'''

Alpha_lineages = ['B.1.1.7', 'B.1.1.28',  'Q.']
Beta_lineages = ['B.1.351', 'B.1.393', 'B.1.416', 'B.1.425', 'B.1.428.2']
Gamma_lineages = ['P.1']
Delta_lineages = ['B.1.617.2', 'AY.']
Epsilon_lineages = ['B.1.427', 'B.1.429', 'B.1.466.2', 'B.1.470', 'B.1.517']
Eta_lineages = ['B.1.525']
Zeta_lineages = ['P.2']
Theta_lineages = ['P.3', 'P.7', 'B.1.1.33', 'B.1.1.97', 'B.1.1.117', 'B.1.1.161', 'B.1.1.216',
                  'B.1.1.222', 'B.1.1.263', 'R.1', 'B.1.1.317', 'B.1.1.376', 'B.1.1.409', 'B.1.1.486']
Iota_lineages = ['B.1.526', 'B.1.558', 'B.1.567',
                 'B.1.575', 'B.1.577', 'B.1.595', 'B.1.609']
Kappa_lineages = ['B.1.617.1']
Lambda_lineages = ['C.37']
Mu_lineages = ['B.1.621', 'B.1.621.1']
BA1_lineages = ['BA.1', 'B.1.1.529', 'BC']
BA2_lineages = ['BA.2', 'BG', 'BH', 'BJ', 'BL', 'BM', 'BN', 'BP', 'BR', 'BS', 'BY',
                'CA', 'CB', 'CH', 'CJ', 'CM', 'CV', 'DD', 'DV', 'FR']
BA3_lineages = ['BA.3']
BA4_lineages = ['BA.4', 'CS', 'DC']
BA5_lineages = ['BA.5', 'BE', 'BF', 'BK', 'BQ', 'BT', 'BU', 'BV', 'BW', 'BZ', 'CD', 'CE',
                'CF', 'CG', 'CK', 'CL', 'CN', 'CP', 'CQ', 'CR', 'CT', 'CU', 'CW', 'CY', 'CZ',
                'DA', 'DB', 'DE', 'DF', 'DG', 'DH', 'DJ', 'DK', 'DL', 'DR', 'DT', 'EY']


XBB_lineages = ['XBB', 'EG', 'EK', 'EL', 'EM', 'EU', 'FD', 'FE', 'FG', 'FG',
                'FH', 'FL', 'FP', 'FT', 'FU', 'FW', 'FY', 'FZ', 'GA', 'GB',
                'GC', 'GD', 'GE', 'GF', 'GG', 'GJ', 'GK']

# Numeric_lineages can house ...numeric lineages. Note: 1.617.3 is a placeholder. Does not exist in AZ dataset.
Numeric_lineages = ['1.617.3', '1.1.318', 'AZ', 'HZ']
Recombinant_lineages = ['XA', 'XB', 'XD', 'XE', 'XF', 'XG', 'XH', 'XJ', 'XK', 'XL',
                        'XM', 'XN', 'XP', 'XQ', 'XR', 'XS', 'XT', 'XU', 'XV', 'XW',
                        'XY', 'XZ', 'XAA', 'XAB', 'XAC', 'XAD', 'XAE', 'XAF', 'XAG',
                        'XAH', 'XAJ', 'XAK', 'XAL', 'XAM', 'XAN', 'XAP', 'XAQ', 'XAR',
                        'XAS', 'XAT', 'XAU', 'XAV', 'XAW', 'XAY', 'XAZ', 'XBA', 'XBC',
                        'XBD', 'XBE', 'XBF', 'XBG', 'XBH', 'XBJ', 'XBK', 'XBL', 'XBM',
                        'XBN', 'XBP', 'XBQ', 'XBR', 'XBS', 'XBT', 'XBU', 'XBV', 'XBW',
                        'XBY', 'XBZ', 'XCA', 'XCB', 'XCC', 'XCF']

'''
FFF1
Adjust the two lists below to reflect accurate sum of total variatn lists above.
All_lineages = list of each individual lineage
Note: Common errors from this section arise due to missing commas or + signs.
'''

All_lineages = Alpha_lineages + Beta_lineages + Gamma_lineages + Delta_lineages + Epsilon_lineages + Eta_lineages \
    + Zeta_lineages + Theta_lineages + Iota_lineages + Kappa_lineages + Lambda_lineages + Mu_lineages + BA1_lineages \
    + BA2_lineages + BA3_lineages + BA4_lineages + BA5_lineages + \
    Numeric_lineages + Recombinant_lineages + XBB_lineages  # U002

All_lineages_list = [Alpha_lineages, Beta_lineages, Gamma_lineages, Delta_lineages, Epsilon_lineages, Eta_lineages,
                     Zeta_lineages, Theta_lineages, Iota_lineages, Kappa_lineages, Lambda_lineages, Mu_lineages,
                     BA1_lineages, BA2_lineages, BA3_lineages, BA4_lineages, BA5_lineages, Numeric_lineages,
                     Recombinant_lineages, XBB_lineages]
'''
GGG1
The three lists below should contain lineages respective to their designation.
Freely add or remove lineages as needed. Adjust repective dictionaries if name change is necessary.
'''
og_VOC_list = ['B.1.1.7', 'B.1.351', 'P.1', 'B.1.617.2', 'B.1.427', 'B.1.525', 'P.2', 'P.3',
               'B.1.526', 'B.1.617.1', 'C.37', 'B.1.621', 'BA.1', 'BA.2', 'BA.3', 'BA.4', 'BA.5', 'XBB']

VOC_list = ['XBB.1.5', 'BQ.1.1', 'BQ.1', 'BA.5', 'XBB', 'BN.1', 'BF.7', 'BA.2.75',
            'BA.5.2.6', 'BA.4.6', 'BA.2', 'BF.11', 'BA.2.75.2', 'BA.4', 'BA.1.1',
            'B.1.1.529', 'BA.2.12.1']

VOI_list = []

VBM_list = ['B.1.1.7',
            'B.1.351',
            'P.1',
            'B.1.617.2',
            'B.1.427',
            'B.1.525'
            'B.1.526',
            'B.1.617.1'
            'P.2'
            'B.1.621'
            ]

'''
VOC_dict is used to rename lineage columns
The dictionaries below should contain all lineages that require name change
Check here first for inconsistencies in nameing.
New entries can be added to the end of the dictionary, old entries can be removded to reduce clutter, but is not necessary.
Overlap in naming designation should result in the last instance taking priority (untested).
'''
VOC_dict = {'B.1.1.7': 'Alpha (B.1.1.7)', 'B.1.351': 'Beta (B.1.351)', 'P.1': 'Gamma (P.1)', 'B.1.617.2': 'Delta (B.1.617.2)',
            'B.1.427': 'Epsilon (B.1.427/429)', 'B.1.525': 'Eta (B.1.525)', 'B.1.526': 'Iota (B.1.526)', 'P.2': 'Zeta (P.2)',
            'B.1.617.1': 'Kappa (B.1.617.1)', 'C.37': 'Lambda (C.37)',  'B.1.621': 'Mu (B.1.621)', 'BA.1': 'Omicron (BA.1/B.1.1.529)', 'BA.2': 'Omicron (BA.2)',
            'BA.2.12.1': 'Omicron (BA.2.12.1)', 'BA.2.75': 'Omicron (BA.2.75)', 'BA.3': 'Omicron (BA.3)', 'BA.4': 'Omicron (BA.4)',
            'BA.5': 'Omicron (BA.5)', 'BQ.1.1': 'Omicron_BA.5.3.1.1.1.1.1.1 (BQ.1.1)', 'BQ.1': 'Omicron_BA.5.3.1.1.1.1.1 (BQ.1)',
            'XBB': 'Recombinant_BA.2 (XBB)', 'BN.1': 'Omicron_BA.2.75.51 (BN.1)', 'BF.7': 'Omicron_BA.5.2.11 (BF.7)',
            'BA.5.2.6': 'Omicron (BA.5.2.6)', 'BA.4.6': 'Omicron (BA.4.6)', 'BF.11': 'Omicron_BA.5.2.1.11 (BF.11)',
            'BA.2.75.2': 'Omicron (BA.2.75.2)', 'BA.1.1': 'Omicron (BA.1.1)', 'B.1.1.529': 'Omicron (B.1.1.529)', 'XBB': 'Recombinant (XBB)'}

# VOI_dict is only used to rename VOI columns
VOI_dict = {'B.1.617.2': 'Delta (B.1.617.2)'}

# VBM_dict1 is only used to rename VBM1 columns. VBM_dict1 is just a more complete version of VBM_dict2. VBM_dict2 still gets the job done.
VBM_dict1 = {'Alpha': [Alpha_lineages], 'Beta': [Beta_lineages], 'Gamma': [Gamma_lineages], 'Delta': [Delta_lineages],
             'Epsilon': ['B.1.427', 'B.1.429'], 'Eta': ['B.1.525'], 'Iota': ['B.1.526'], 'Kappa': ['B.1.617.1'],
             'Numberic': ['1.617.3'], 'Mu': ['B.1.621', 'B.1.621.1'], 'Zeta': ['P.2']}

# VBM_dict2 is only used to rename VBM2 columns
VBM_dict2 = {'B.1.1.7': 'Alpha (B.1.1.7)', 'B.1.617.2': 'Delta (B.1.617.2)', 'P.1': 'Gamma (P.1)', 'B.1.351': 'Beta (B.1.351)',
             'B.1.427': 'Epsilon (B.1.427/429)', 'B.1.525': 'Eta (B.1.525)', 'B.1.526': 'Iota (B.1.526)', 'P.2': 'Zeta (P.2)',
             'B.1.617.1': 'Kappa (B.1.617.1)', 'B.1.621': 'Mu (B.1.621)', 'BA.1': 'Omicron (BA.1/B.1.1.529)', 'BA.2': 'Omicron (BA.2)',
             'BA.3': 'Omicron (BA.3)', 'BA.4': 'Omicron (BA.4)', 'BA.5': 'Omicron (BA.5)'}


'''
HHH1
'''
primary_mutation_list = ['Spike_R346T', 'Spike_K444T', 'Spike_N460K',
                         'Spike_F486S', 'Spike_K444M', 'Spike_L452R', 'Spike_F486V']
# secondary_mutation_list = []
# lineage_list = ['BQ.1.1', 'BA.2.75.2', 'BA.5.2.7', 'BR.1']
# secondary_lineage_list = []
# res = []

'''
III1
Color number to change color of output terminal text.
'''
color_number = 92

'''
JJJJ
Select VOC output: 1) VOC_list, og_VOC_list
'''
VOC_output = og_VOC_list

#######################################################################################################################
# Function exicutuion
#######################################################################################################################

'''
File formats:
infile1 = 'CombinedGISAID_hcov-19_2022_12_20_18.tsv'
infile2 = 'gisaid_hcov-19_2022_12_27_22.tsv'
outfile1 = 'CombinedGISAID_hcov-19_2022_12_27_22.tsv'
Note: these dont get updated automatically.
'''

gisaid_update(infile1, infile2).to_csv(
    outfile1, mode='w', sep='\t', header=True)

# df contains the columns selected from target csv, using column list t1.
sorted_t1_data_df = gisaid_data_sort_t1(outfile1, t1_columns)
# sorted_t1_data_df.to_csv(outfile4)  # outfile4 = f'{outdate}_GISAID_t1_output.csv'

# df contains the output of t2 sort using the sorted_t1_data_df, columns from t2a list,
# the string for the column being counted, and the string for the sorting method !!case sensitive!! 'Week' or 'Day'.
sorted_t2a_data_df = gisaid_data_sort_t2(
    sorted_t1_data_df, t2a_columns, 'Lineage', 'Week')  # EEEE
# outfile2 = f'{outdate}_lineage_by_week_report.csv'
sorted_t2a_data_df.to_csv(outfile2)

'''
The following dfs after df01 house the lineage_count function outputs for each lineage list within the funciton.
df output should be a single column named after 'designation'_lineages[0] that contains the count for all lineages within
the specific lineage list.
eg. df02 contains the
'''
df = sorted_t2a_data_df.copy()
df01 = pd.DataFrame(df.index).set_index('Week')
df02 = pd.DataFrame(lineage_count(Alpha_lineages, sorted_t2a_data_df))  # DDDD
df03 = pd.DataFrame(lineage_count(Beta_lineages, sorted_t2a_data_df))
df04 = pd.DataFrame(lineage_count(Gamma_lineages, sorted_t2a_data_df))
df05 = pd.DataFrame(lineage_count(Delta_lineages, sorted_t2a_data_df))
df06 = pd.DataFrame(lineage_count(Epsilon_lineages, sorted_t2a_data_df))
df07 = pd.DataFrame(lineage_count(Zeta_lineages, sorted_t2a_data_df))
df08 = pd.DataFrame(lineage_count(Theta_lineages, sorted_t2a_data_df))
df09 = pd.DataFrame(lineage_count(Iota_lineages, sorted_t2a_data_df))
df10 = pd.DataFrame(lineage_count(Lambda_lineages, sorted_t2a_data_df))
df11 = pd.DataFrame(lineage_count(Mu_lineages, sorted_t2a_data_df))
df12 = pd.DataFrame(lineage_count(BA1_lineages, sorted_t2a_data_df))
df13 = pd.DataFrame(lineage_count(BA2_lineages, sorted_t2a_data_df))
df14 = pd.DataFrame(lineage_count(BA3_lineages, sorted_t2a_data_df))
df15 = pd.DataFrame(lineage_count(BA4_lineages, sorted_t2a_data_df))
df16 = pd.DataFrame(lineage_count(BA5_lineages, sorted_t2a_data_df))
df17 = pd.DataFrame(lineage_count(Numeric_lineages, sorted_t2a_data_df))
df18 = pd.DataFrame(lineage_count(Recombinant_lineages, sorted_t2a_data_df))
df_final = pd.concat([df02, df03, df04, df05, df06, df07, df08,
                      df09, df10, df11, df12, df13, df14, df15, df16, df17, df18], axis=1)

# df02.to_csv('df02_test.csv')
# df_final.to_csv('df_final_test.csv')

#######################################################################################################################
# Summary Tables
#######################################################################################################################

df_final_a = df_final.copy()
for i in VOC_output:
    try:
        df01[i] = df[i]
    except Exception:
        pass
df01.update(df_final_a)
dfa = column_drop(df, All_lineages)

# Dropping samples that have blank lineage designations
# print("Dropping samples without lineage names. This contains the following information: ")
# df = df.drop(columns=['NA'])

# Count all remaining dataframes and consider them "Others"
df01['Other'] = 0
print("Creating other variant column, this contains the following lineages: ")
col_list = []
for col in dfa:
    col_list.append(col)
    df01['Other'] += dfa[col]
print(col_list)
df01.rename(columns=dict(VOC_dict), inplace=True)
df01.to_csv('VariantTable.csv', index=True)

sorted_t2b_data_df = gisaid_data_sort_t2(
    sorted_t1_data_df, t2b_columns, 'Lineage', 'Ouptput_date_format')  # CCCC
tb2_date_summary_df = gisaid_date_summary_count(sorted_t2b_data_df)  # EEEE
tb2_date_summary_df.to_csv(outfile3)

dfb = df.copy()
df_final_b = df_final.copy()
dfb.update(df_final_b)
summary = gisaid_date_summary_count(dfb)

print(color_me('Printing: VOC Summary', color_number))
print(table_summary(summary, VOC_output, VOC_dict))  # GGGG
table_summary(summary, VOC_output, VOC_dict).to_csv(
    'VOCTableSummary.csv', index=True)
print(color_me('Printing: VOI_Summary', color_number))
print(table_summary(summary, VOI_list, VOI_dict))
table_summary(summary, VOI_list, VOI_dict).to_csv(
    'VOITableSummary.csv', index=True)
print(color_me('Printing: VBM_summary', color_number))
print(table_summary(summary, VBM_list, VBM_dict2))
table_summary(summary, VBM_list, VBM_dict2).to_csv(
    'VBMTableSummary.csv', index=True)

mutation_counter(outfile1, primary_mutation_list)  # HHHH

sys.stdout.close()
sys.stdout = stdoutOrigin

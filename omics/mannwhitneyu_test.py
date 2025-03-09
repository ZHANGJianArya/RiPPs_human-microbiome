import pandas as pd
import scipy.stats as stats

def check_normality(data):
    """
    Check the normality of the data using the Shapiro-Wilk test.

    Parameters:
    data (array-like): The data to test for normality.

    Returns:
    float: The test statistic for the Shapiro-Wilk test.
    float: The p-value for the Shapiro-Wilk test.
    """
    test_stat_normality, p_value_normality = stats.shapiro(data)
    return test_stat_normality, p_value_normality

# Read the input data file
df = pd.read_csv("RiPP_TPM_IBD_MT_cluster90_p12_stat.tsv", encoding="UTF-8", engine="python", sep="\t", index_col=0)

# Get a list of all columns in the dataframe without the Group column
column_list = [x for x in df.columns if x != 'Group']

# Create empty dictionaries to store test results
tuc_test_results = {}
tcd_test_results = {}

# Loop over column_list and perform statistical tests
for column in column_list:
    group_ctr = df.where(df.Group == "nonIBD").dropna()[column]
    group_uc = df.where(df.Group == "UC").dropna()[column]

    test_stat_normality_ctr, p_value_normality_ctr = stats.shapiro(group_ctr)
    test_stat_normality_uc, p_value_normality_uc = stats.shapiro(group_uc)

    if p_value_normality_ctr > 0.05 and p_value_normality_uc > 0.05:
        tuc_test_results[column] = stats.ttest_ind(group_ctr, group_uc)
    else:
        tuc_test_results[column] = stats.mannwhitneyu(x=group_ctr, y=group_uc, alternative='two-sided')

# Convert results to DataFrame and save to CSV
results_uc_df = pd.DataFrame.from_dict(tuc_test_results, orient='Index')
results_uc_df.columns = ['statistic', 'pvalue']
results_uc_df.to_csv("RiPP_in_IBD_uc_MT_sta_tpm.tsv", sep="\t", index=True)

for column in column_list:
    group_ctr = df.where(df.Group == "nonIBD").dropna()[column]
    group_cd = df.where(df.Group == "CD").dropna()[column]

    test_stat_normality_ctr, p_value_normality_cd = stats.shapiro(group_ctr)
    test_stat_normality_cd, p_value_normality_cd = stats.shapiro(group_cd)

    if p_value_normality_ctr > 0.05 and p_value_normality_cd > 0.05:
        tcd_test_results[column] = stats.ttest_ind(group_ctr, group_cd)
    else:
        tcd_test_results[column] = stats.mannwhitneyu(x=group_ctr, y=group_cd, alternative='two-sided')

# Convert results to DataFrame and save to CSV
results_cd_df = pd.DataFrame.from_dict(tcd_test_results, orient='Index')
results_cd_df.columns = ['statistic', 'pvalue']
results_cd_df.to_csv("RiPP_in_IBD_cd_MT_sta_tpm.tsv", sep="\t", index=True)

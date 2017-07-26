# A thorough analysis of a specific dataset
# For debugging and method comparison
dataset = "GSE18608;tissue: peripheral blood;GPL570" # Two time points, no sex, age
dataset = "GSE28998;tissue: biceps;GPL570" # Two time points with sex
dataset = "GSE59088;tissue: Vastus lateralis;GPL6244" # More than two time points
# A dataset with sex and age
dataset_samples = names(which(analysis_dataset_ids==dataset))
# cut by samples that have time points and are not repeats
dataset_samples = dataset_samples[metadata[dataset_samples,time_col] != ""]
dataset_samples = dataset_samples[metadata[dataset_samples,"Replicate.info"] == ""]
curr_platforms = as.character(metadata[dataset_samples,"platform_id"])
platform = curr_platforms[1]

# Filter the current sample set:
curr_samples = dataset_samples
curr_samples = curr_samples[sample2replicate_info[curr_samples] == ""]
curr_samples = curr_samples[sample2subject[curr_samples]!=""]
curr_times = as.character(metadata[curr_samples,time_col])
curr_samples = curr_samples[curr_times !=""]

# use frma
data_matrix = get_data_matrix_from_matrix_lists(CEL_frma_profiles,curr_samples)
curr_samples = curr_samples[order(metadata[curr_samples,time_col])]
curr_times = as.numeric(as.character(metadata[curr_samples,time_col]))
curr_subjects = factor(sample2subject[curr_samples])
curr_ages = factor(sample2age[curr_samples],ordered=T)
curr_sex = factor(sample2sex[curr_samples])

# transform to genes
genes_data_matrix_obj = transform_matrix_into_genes(data_matrix,gpl = platform,gpl_mappings_entrez2probes)
genes_data_matrix = genes_data_matrix_obj$entrez_mat

# fold changes
gene_fold_changes = get_fold_changes_vs_baseline(genes_data_matrix,curr_subjects,as.numeric(curr_times))

# Look at the top fchange genes
fc_thr = 1
gene_fc_count = apply(gene_fold_changes>fc_thr,1,sum)
gene_fc_count_abs = apply(abs(gene_fold_changes)>fc_thr,1,sum)
gene_fc_count_minus = apply(gene_fold_changes< -fc_thr,1,sum)
plot(gene_fc_count,gene_fc_count_minus)

# select genes
genes = gene_fc_count >= max(sort(gene_fc_count,decreasing=T)[1000],2)
genes_data_matrix = genes_data_matrix[genes,]
print(dim(genes_data_matrix))
ord = order(curr_subjects)
gene_profile = genes_data_matrix[1,]
gene_fc_profile = gene_fold_changes[rownames(genes_data_matrix)[1],]
plot(gene_profile[ord],pch = as.numeric(factor(curr_subjects[ord])))

run_paired_test(gene_profile,curr_subjects,curr_times)$p.value
run_paired_test(gene_profile,curr_subjects,curr_times,func=wilcox.test)$p.value
run_anova_on_time_and_subject(gene_profile,curr_subjects,curr_times)[[2]]
library(lme4)
run_lmer4_anova(gene_profile,curr_subjects,curr_times)
# with sex and age
run_lmer4_anova(gene_profile,curr_subjects,curr_times,sexv=curr_sex,agev = curr_ages)
run_anova_on_time_age_sex_and_subject(gene_profile,curr_subjects,curr_times,sexv=curr_sex,agev=curr_ages)[[2]]

# run all three methods and compare
curr_paired_tests = apply(genes_data_matrix,1,run_anova_on_time_age_sex_and_subject,
  subjs=curr_subjects,timev=curr_times,sexv=NULL,agev=NULL)
gene_paired_test_pvals1 = sapply(curr_paired_tests,function(x)x$pval)
curr_paired_tests2 = apply(genes_data_matrix,1,run_paired_test,subjs=curr_subjects,timev=curr_times)
gene_paired_test_pvals2 = sapply(curr_paired_tests2,function(x)x$p.value)
gene_paired_tests_lme4_pvals = apply(genes_data_matrix,1,run_lmer4_anova,
                   subjs=curr_subjects,timev=curr_times,sexv=NULL,agev=NULL)
curr_paired_tests3 = apply(genes_data_matrix,1,run_anova_on_time_age_sex_and_subject,
                           subjs=curr_subjects,timev=curr_times,sexv=curr_sex)
gene_paired_test_pvals3 = sapply(curr_paired_tests3,function(x)x$pval)
gene_paired_tests_lme4_pvals = apply(genes_data_matrix,1,run_lmer4_anova,
                                     subjs=curr_subjects,timev=curr_times,sexv=curr_sex,agev=NULL)

plot(gene_paired_test_pvals1,gene_paired_test_pvals2)
plot(gene_paired_test_pvals1,gene_paired_test_pvals3);abline(0,1)
plot(gene_paired_test_pvals1,gene_paired_tests_lme4_pvals);abline(0,1)
plot(gene_paired_test_pvals3,gene_paired_tests_lme4_pvals);abline(0,1)
plot(sort(p.adjust(gene_paired_tests_lme4_pvals)))

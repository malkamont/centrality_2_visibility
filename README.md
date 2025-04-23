# centrality_2_visibility
Replication data + code for "Beyond Formal Power: How Central Roles in Political Networks are Related to Media Visibility".

1. Run extract_visibility.R [scripts] to retrieve the media appearance counts of organisations from a dna object. NOTE: due to data protection, we provide only hard copies of these counts.
2. Run reconstruct_collaboration.R in the arr environment [yamls] to reconstruct the two collaboration networks. NOTE: due to a version/dependency mismatch between ergm and Bergm packages, fitted objects [hard_copies] should be used for exact reproducibility.
3. Run data_processing.py in the gt environment as well as climate_progressiveness.R and data_processing.R in the arr environment to process the data. NOTE that due to data protection, we provide the final joint data only as a hard copy.
4. Run data_describe.R in the arr environment to collect descriptive statistics on the final joint data.
5. Run model_fit.R in the meccess environment to fit the Bayesian regression models. Access to a high-performance computing cluster speeds up the process.
6. Run model_summary.R in the meccess environment to summarise and illustrate the results.

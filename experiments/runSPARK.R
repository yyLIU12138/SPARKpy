library('SPARK')

save_dir <- './results/Rep11_MOB/'
# i <- 0

rawcount <- read.csv(paste0(save_dir, "data_for_SPARK/expression_matrix.csv"), row.names = 1)
# rawcount <- read.csv(paste0(save_dir, "data_for_SPARK/expression_matrix_", i, ".csv"), row.names = 1)

colnames(rawcount) <- gsub("^X", "", colnames(rawcount))

info <- read.csv(paste0(save_dir, "data_for_SPARK/spatial_coordinates.csv"), row.names = 1)

## filter genes and cells/spots and 
spark <- CreateSPARKObject(counts=rawcount, 
                           location=info[,1:2],
                           percentage = 0.1, 
                           min_total_counts = 10)

## total counts for each cell/spot
spark@lib_size <- apply(spark@counts, 2, sum)

## Estimating Parameter Under Null
spark <- spark.vc(spark, 
                  covariates = NULL, 
                  lib_size = spark@lib_size, 
                  num_core = 1,
                  verbose = F)

## Calculating pval
spark <- spark.test(spark, 
                    check_positive = T, 
                    verbose = F)

# write.csv(spark@res_mtest[, c("GSP1", "COS1", "GSP2", "COS2", 
#                               "GSP3", "COS3", "GSP4", "COS4", 
#                               "GSP5", "COS5", "combined_pvalue", "adjusted_pvalue")], 
#           file = paste0(save_dir, "SPARK_R_p_values_permute_", i, ".csv"), 
#           row.names = TRUE)

write.csv(spark@res_mtest[, c("GSP1", "COS1", "GSP2", "COS2",
                              "GSP3", "COS3", "GSP4", "COS4",
                              "GSP5", "COS5", "combined_pvalue", "adjusted_pvalue")],
          file = paste0(save_dir, "SPARK_R_p_values.csv"),
          row.names = TRUE)



# download data

suppressPackageStartupMessages({
  library(HDCytoData)
})

AML1 <- Weber_AML_sim_main_0.1pc_SE()
AML2 <- Weber_AML_sim_main_1pc_SE()
AML3 <- Weber_AML_sim_main_5pc_SE()

saveRDS(object = AML1, file = "data/Weber_AML_sim_main_0.1pc_SE.rds")
saveRDS(object = AML2, file = "data/Weber_AML_sim_main_1pc_SE.rds")
saveRDS(object = AML3, file = "data/Weber_AML_sim_main_5pc_SE.rds")


BCR1 <- Weber_BCR_XL_sim_less_distinct_less_75pc_SE()
BCR2 <- Weber_BCR_XL_sim_less_distinct_less_50pc_SE()
BCR3 <- Weber_BCR_XL_sim_main_SE()

saveRDS(object = BCR1, file = "data/Weber_BCR_XL_sim_less_distinct_less_75pc_SE.rds")
saveRDS(object = BCR2, file = "data/Weber_BCR_XL_sim_less_distinct_less_50pc_SE.rds")
saveRDS(object = BCR3, file = "data/Weber_BCR_XL_sim_main_SE.rds")

library(lme4)        
library(lmerTest)    
library(boot)
modalities <- c("both", "face", "word")
sources <- c("subavg", "subspec")

slurm_array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
modality_idx <- (slurm_array_id - 1) %/% length(sources) + 1
source_idx   <- (slurm_array_id - 1) %% length(sources) + 1
modalities <- c(modalities[modality_idx])
sources <- c(sources[source_idx])

df_all <- read.csv("./outputs/inoutscan_beh_allsubs_wMaxCategoryDistance.csv")
for (modality in modalities) {
  for (source in sources) {
    print("--------------------------------")
    print(source)
    print(modality)
    if (modality == "both") {
      df <- df_all[df_all$source == source, ]
    } else {
      df <- df_all[df_all$source == source & df_all$modality == modality, ]
    }
    #drop nan in outscan_rating_change and outscan_consistency
    df <- df[!is.na(df$outscan_rating_change) & !is.na(df$outscan_consistency), ]
    df$sub_id <- as.factor(df$sub_id)
    df$outscan_consistency <- as.numeric(df$outscan_consistency)
    #df$outscan_rating_change_c <- scale(df$outscan_rating_change, scale = FALSE)#mean center
    #df$outscan_rating_change_z <- scale(df$outscan_rating_change) #z

    m <- glmer(outscan_consistency ~ outscan_rating_change + 
                  (1 + outscan_rating_change | sub_id),
                data = df, family = binomial)
    print(summary(m))

    m2 <- glmer(outscan_consistency ~ outscan_rating_change + 
                                          I(outscan_rating_change^2) +
                                          (1 + outscan_rating_change | sub_id),
                  data = df, family = binomial)
    print(anova(m, m2))
    print(summary(m2))


    compute_estimate <- function(fit) {fixef(fit)["outscan_rating_change"]}

    set.seed(123)
    boot_results <- bootMer(m,FUN = compute_estimate,nsim = 10000,type = "parametric")
    ci <- quantile(boot_results$t, c(0.025, 0.975))
    print('CI for m:')
    print(ci)
    boot_results <- bootMer(m2,FUN = compute_estimate,nsim = 10000,type = "parametric")
    ci <- quantile(boot_results$t, c(0.025, 0.975))
    print('CI for m2:')
    print(ci)
  }
}

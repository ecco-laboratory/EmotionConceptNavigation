library(lme4)        
library(lmerTest)    
library(boot)
library(dplyr)
project_dir = '/home/data/eccolab/MNS/';


region_names <- c("OFC2016ConstantinescuR5","HC", "aHC","pHC")
sources <- c("Subavg", "Subspec")
distance_types <- c("va", "judged", "catCosine", "catMax",'catDominantAvg','judgedSigned')
modality_to_use = c('both_modalities');


slurm_array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
type_idx     <- ((slurm_array_id - 1) %% length(distance_types)) + 1
source_idx   <- (((slurm_array_id - 1) %/% length(distance_types)) %% length(sources)) + 1
region_idx   <- (((slurm_array_id - 1) %/% (length(distance_types) * length(sources))) %% length(region_names)) + 1
modality_idx <- (((slurm_array_id - 1) %/% (length(distance_types) * length(sources) * length(region_names))) %% length(modality_to_use)) + 1
regions <- c(region_names[region_idx])
sources <- c(sources[source_idx])
distance_types <- c(distance_types[type_idx])
modality_to_use <- c(modality_to_use[modality_idx])
print('multivariate encoding')
print('controlling for modality')
cat("SLURM_ARRAY_TASK_ID:", slurm_array_id, "\n")
cat("Modality:", modality_to_use, "\n")
cat("Region:", regions, "\n")
cat("Source:", sources, "\n")
cat("Distance type:", distance_types, "\n")


for (modality in modality_to_use) {
  for (region in regions) {
    for (source in sources) {
      for (distance_type in distance_types) {

        df <- read.csv(paste0(
              project_dir,
              'outputs/singleTrialBetaAnalysis/noICA/incl_all_subs_trials/',
              'distanceCentered_multivariate_encoding_yhat_ytest/unsmoothed/',
              modality, '/', source, '/', distance_type, '/csv/',
              'beh_Yproj_raw_wDistance_', region, '.csv'
            ))

        df <- df[!is.na(df$Yproj) & !is.na(df$consistency), ]
        df$test_subject <- as.factor(df$test_subject)
        df$consistency <- as.numeric(df$consistency)
        df$Yproj <- df$Yproj - mean(df$Yproj);


        df <- df %>%
          group_by(test_subject) %>%
          mutate(modality = c(rep("word", 124),
                              rep("face", 84))) %>%
          ungroup()

        df$modality <- as.factor(df$modality)

        
        #this is the model to test if the effect of Yproj's relationship with consistency is different for word and face
        m <- glmer(consistency ~ Yproj + modality + Yproj * modality + (1 + Yproj + modality | test_subject),
                    data = df, family = binomial)
        print(summary(m))

        m2 <- glmer(consistency ~ Yproj + (1 + Yproj | test_subject),
                   data = df, family = binomial)
        print(anova(m, m2))
        print(summary(m2))

        compute_estimate <- function(fit) {fixef(fit)}

        set.seed(123)
        boot_results <- bootMer(m,FUN = compute_estimate,nsim = 10000,type = "parametric")
        ci <- apply(boot_results$t, 2, quantile, probs = c(0.025, 0.975))
        print('CI for m:')
        print(ci)
      }
    }
  }
}


library(data.table)
library(ggplot2)
# 0.00089325591782
readLines(
  "AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210.pro.cluster.txt"
)[44]
sip0 <-
  fread(
    "AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210.pro.cluster.txt",
    skip = 81
  )
# 2241
nrow(sip0)
# 0.1162298
mean(sip0$AverageEnrichmentLevel)
# 0.4
median(sip0$AverageEnrichmentLevel)

# 0.0
readLines(
  "AMD_StandardEnrichment_SampleBeta_50Percent15N_Velos_OrbiMS2_Run2_020110.pro.cluster.txt"
)[44]
sip50 <-
  fread(
    "AMD_StandardEnrichment_SampleBeta_50Percent15N_Velos_OrbiMS2_Run2_020110.pro.cluster.txt",
    skip = 81
  )
# 853
nrow(sip50)
# 46.89325
mean(sip50$AverageEnrichmentLevel)
# 45
median(sip50$AverageEnrichmentLevel)

# 0.000478240076518
readLines(
  "AMD_StandardEnrichment_SampleAlpha_98Percent15N_Velos_OrbiMS2_Run2_013110.pro.cluster.txt"
)[44]
sip98 <-
  fread(
    "AMD_StandardEnrichment_SampleAlpha_98Percent15N_Velos_OrbiMS2_Run2_013110.pro.cluster.txt",
    skip = 81
  )
# 2092
nrow(sip98)
# 99.63204
mean(sip98$AverageEnrichmentLevel)
# 100
median(sip98$AverageEnrichmentLevel)

x <-
  data.frame(
    Abundance = c(
      sip0$AverageEnrichmentLevel,
      sip50$AverageEnrichmentLevel,
      sip98$AverageEnrichmentLevel
    ),
    SIP = c(rep(0, nrow(sip0)), rep(50, nrow(sip50)), rep(98, nrow(sip98)))
  )
x$SIP <- as.factor(x$SIP)
p <-
  ggplot(data = x,
         mapping = aes(x = Abundance, fill = SIP)) +
  geom_histogram(binwidth = 1,
                 color = I("black")) +
  xlab("N15 abundance (%)") +
  ylab("Protein Count") +
  theme(text = element_text(size = 15))
p + ggsave("SIPresult.pdf", width = 8, height = 6)

library(data.table)
library(ggplot2)
# 0.0
readLines(
  "AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210.pro.cluster.txt"
)[44]
sip0 <-
  fread(
    "AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210.pro.cluster.txt",
    skip = 81
  )
# 1887
nrow(sip0)
# 5.052384
mean(sip0$AverageEnrichmentLevel)
# 5
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
# 836
nrow(sip50)
# 48.70662
mean(sip50$AverageEnrichmentLevel)
# 48.33333
median(sip50$AverageEnrichmentLevel)

# 0.000
readLines(
  "AMD_StandardEnrichment_SampleAlpha_98Percent15N_Velos_OrbiMS2_Run2_013110.pro.cluster.txt"
)[44]
sip98 <-
  fread(
    "AMD_StandardEnrichment_SampleAlpha_98Percent15N_Velos_OrbiMS2_Run2_013110.pro.cluster.txt",
    skip = 81
  )
# 2007
nrow(sip98)
# 94.93589
mean(sip98$AverageEnrichmentLevel)
# 95
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

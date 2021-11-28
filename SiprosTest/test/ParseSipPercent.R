library(data.table)
library(ggplot2)
sip0 <-
  fread(
    "AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210.pro.cluster.txt",
    skip = 79
  )
# 4025
nrow(sip0)
# 0.6467496
mean(sip0$AverageEnrichmentLevel)
# 0.3918919
median(sip0$AverageEnrichmentLevel)
sip50 <-
  fread(
    "AMD_StandardEnrichment_SampleBeta_50Percent15N_Velos_OrbiMS2_Run2_020110.pro.cluster.txt",
    skip = 79
  )
# 46.88082
mean(sip50$AverageEnrichmentLevel)
# 47.25
median(sip50$AverageEnrichmentLevel)
sip98 <-
  fread(
    "AMD_StandardEnrichment_SampleAlpha_98Percent15N_Velos_OrbiMS2_Run2_013110.pro.cluster.txt",
    skip = 79
  )
# 95.29632
mean(sip98$AverageEnrichmentLevel)
# 96.83935
median(sip98$AverageEnrichmentLevel)

x <-
  data.frame(
    Abundance = c(sip0$AverageEnrichmentLevel,
                  sip50$AverageEnrichmentLevel,sip98$AverageEnrichmentLevel),
    SIP = c(rep(0, nrow(sip0)), rep(50, nrow(sip50)),rep(98, nrow(sip98)))
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


library(ggthemes)
library(tidyverse)
load('results/MSG.master_meta.Rdata')

meta = meta %>% mutate(NM = factor(ifelse(grepl("NM", HB), "NM", "M"), levels = c("NM", "M")), HB = gsub("-NM", "", HB))

meta %>% filter(TimePoint != "Recurrence 3" & TimePoint != "Recurrence 2") %>%
  ggplot(aes(x = Pt)) + 
  geom_bar(aes(fill = HB, alpha=NM)) + 
  facet_wrap( ~ TimePoint, ncol = 1) + 
  scale_fill_ptol() +
  theme_fivethirtyeight(base_size = 18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(fill = "", alpha = "", y = "Count", x = "Patient")

meta %>% filter(TimePoint != "Recurrence 3" & TimePoint != "Recurrence 2") %>%
  ggplot(aes(x = Pt)) + 
  geom_bar(aes(fill = LGm, alpha=LGm_prob)) + 
  facet_wrap( ~ TimePoint, ncol = 1) + 
  scale_fill_ptol() +
  theme_fivethirtyeight(base_size = 18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(fill = "", alpha = "", y = "Count", x = "Patient")

meta %>% filter(Dataset == "VUmc") %>%
  ggplot(aes(x = Pt)) + 
  geom_bar(aes(fill = HB)) + 
  facet_wrap( ~ TimePoint, ncol = 1) + 
  scale_fill_ptol() +
  theme_fivethirtyeight(base_size = 18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(fill = "", y = "Count", x = "Patient")

meta %>% filter(Dataset == "VUmc") %>%
  ggplot(aes(x = Pt)) + 
  geom_bar(aes(fill = RF)) + 
  facet_wrap( ~ TimePoint, ncol = 1) +
  scale_fill_ptol() +
  theme_fivethirtyeight(base_size = 18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(fill = "", y = "Count", x = "Patient")

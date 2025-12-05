library(dplyr)
df_filtered <- output[output$prevalence >= 0.05, ]
df_increasing <- df_filtered %>%
  group_by(classified_label) %>%
  arrange(date) %>%
  mutate(increasing = gr > lag(gr)) %>%
  filter(increasing == TRUE)
# 
# df_early_growth <- df_increasing %>%
#   group_by(classified_label) %>%
#   slice(1) %>%
#   ungroup()

write.csv(df_early_growth, "C://Users//young//OneDrive//Desktop//UKHSA_COVID19_Variants//early_growth_uk.csv")

df_increasing$categorized_prevalence <- cut(df_increasing$prevalence,
                                 breaks = c(-Inf, 0.05,0.06, 0.1,Inf),
                                 labels = c("<5%", "5-6%", "7-9%", ">10%"),
                                 #breaks = c(-Inf, 0.05, 0.1,0.15, Inf),
                                 #labels = c("<5%", "5-9%", "10-14%", ">15%"),                                 
                                 right = FALSE) # This means the interval is closed on the right
df_increasing$categorized_prevalence<-factor(df_increasing$categorized_prevalence, levels = c("<5%", "5-6%", "7-9%", ">10%"), 
                                                labels=c("<5%", "5-6%", "7-9%", ">10%"))

mean_growth_rates <- df_increasing %>%
  group_by(classified_label,categorized_prevalence) %>%
  summarise(mean_gr = mean(gr, na.rm = TRUE), 
            mean_gr_lower = min(gr_lower, na.rm = TRUE), 
            mean_gr_upper = max(gr_upper, na.rm = TRUE)
            )%>%
  mutate(country="United Kingdom")

write.csv(mean_growth_rates, "C://Users//young//OneDrive//Desktop//UKHSA_COVID19_Variants//mean_growth_rates1_uk.csv")
print(mean_growth_rates)

sk<-read.csv("C://Users//young//OneDrive//Desktop//UKHSA_COVID19_Variants//mean_growth_rates1_sk.csv")
uk<-read.csv("C://Users//young//OneDrive//Desktop//UKHSA_COVID19_Variants//mean_growth_rates1_uk.csv")
df<-rbind(sk,uk)


all_combinations <- expand.grid(
  classified_label= unique(df$classified_label),
  categorized_prevalence = unique(df$categorized_prevalence),
  country = unique(df$country),
  stringsAsFactors = FALSE
)

# Join with original data to ensure all combinations have entries, filling missing data with NA
complete_df <- left_join(all_combinations, df, by = c("classified_label", "categorized_prevalence", "country"))

complete_df$combined_group <- interaction(complete_df$classified_label, complete_df$categorized_prevalence, sep = " - ")


library(ggplot2)
forest_plot <- ggplot(complete_df, aes(x = mean_gr, y = combined_group , xmin = mean_gr_lower, xmax = mean_gr_upper, color =classified_label)) +
  geom_pointrange() + # Point range for growth rates and confidence intervals
  labs(x = "Growth Rate", y = "Variant - Prevalence") +
  facet_wrap(~country, scales = "free_y", , nrow = 1) + # Faceting by country
  theme_minimal() + # Minimal theme
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + # Rotate x-axis labels for readability
  ggtitle("Growth Rate by Variant and Prevalence across Countries") # Adding a title

forest_plot

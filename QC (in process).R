library(tidyverse)
library(ggsignif) #for significance brackets

#Load raw data files----------------------------------------------------------------
raw_file_names <- "vieRaw_20250221_titration_20250127_plasma_Flux_wo_sky2iso_res.tsv"
raw_data <- read_tsv(raw_file_names)

  

for (metabolite_studied in metabolites_studied) {
  raw_data_QC <- raw_data|>
    filter(str_detect(sample, "QC"))|>  #select QC injections
    filter(metabolite == metabolite_studied)|>
    mutate(log2_area = na_if(log2(area), -Inf))|> #when log2 of 0 interpret NA
    mutate(injection_number = as.character(
      str_extract(sample, "(?<=_)\\d+(?=\\.raw$)")))|>
    mutate(injection_number = as.character(as.integer(injection_number))) |> 
    select(metabolite, isotopologue, log2_area, injection_number)
    #view(raw_data_QC)) checkpoint
  
  #raw_data_QC_filtered <- raw_data_QC |>
   # filter(injection_number %in% c("6", "19")) |>
    #arrange( injection_number, isotopologue)
  
  #paired_data <- raw_data_QC_filtered |>
   # group_by(isotopologue) |>
    #filter(n_distinct(injection_number[!is.na(log2_area)]) == 2) |>
    #ungroup()
  
QC_graph <- ggplot(paired_data, mapping = aes(x = injection_number, y = log2_area))+
  theme_classic() +
  geom_boxplot() 
  geom_signif(
    comparisons = list(c("6", "19")),
    test = "t.test",
    test.args = list(paired = TRUE),
    map_signif_level = TRUE, #adds significance stars
    textsize = 2,
    margin_top = 0.03, #adds space to error bar
    step_increase = 0.1, #adds space between bars
    tip_length = 0.01) +
  labs(  
    title = paste("Area variation across injections"),
    subtitle =paste("Metabolite studied:", each_metabolite),
    x = "injection number",
    y= "log2 area") 
plot(QC_graph) 

  name_QC_graph <- paste("injection date_", each_metabolite, "_QC.png", sep = "")
  ggsave(plot = QC_graph,
         filename = file.path("Plasma", "Analytical control", name_QC_graph), 
         device = "png",
         width = 10,
         height = 6,
         units = "in",
         dpi = 300,
         limitsize = FALSE)
}


  




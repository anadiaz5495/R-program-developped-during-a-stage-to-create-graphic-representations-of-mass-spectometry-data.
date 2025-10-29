library(tidyverse)
library(ggsignif) #for significance brackets
library(mgcv) #for non linear model

#1. IMPORT DOCUMENT AND CREATE DATABASE WITH PERTINENT MODIFICATIONS
  #1.1 Load raw data files----------------------------------------------------------------

print("Please enter file names:")
files <- scan(what= character())
files <- paste0(files, ".tsv")
n_files <- length(files)


raw_data_original <-read_tsv(files, id = "file") #read both files and create a new column for file name
#view(raw_data_original)
  #1.2 Create column with glucose condition, injection number and date from sample names------
raw_data <- raw_data_original |>
  mutate(
    "glucose_condition" = factor(
      str_remove_all(str_extract(sample, "_(5|10|20)_|(5|10|20)mM|(IP-13C-)Glc|(IP-13C-)Lac|(IP-13C-)BHB| (IP-13C-)Ala"), "_|mM|-13C"),
      levels = c("5", "10", "20", "IP-Glc", "IP-Lac", "IP-BHB", "IP-Ala")))|>
  mutate(injection_date = str_extract(sample, "\\d{6}")) |>
  mutate(injection_number = as.numeric(
    str_extract(sample, "(?<=_)\\d+(?=\\.raw$)")))|> #extract number of sample from sample name
  select(file, sample, injection_date, injection_number, metabolite, isotopologue, isotopic_inchi, area, corrected_area, isotopologue_fraction, residuum, mean_enrichment, glucose_condition )
view(raw_data)
metabolites <- unique(raw_data$metabolite)

  #1.3 Select what metabolites to study------------------------------------------------------------
choice_metabolite <- menu( 
  choices = c("All metabolites", "Selection of metabolites"),
  title = "Please choose an option number:")
if (choice_metabolite ==1) {
  metabolites_studied <- metabolites
}

if (choice_metabolite ==2) {
  cat("Please make selection of metabolites to study:")
  metabolites_studied <- scan(what = character(), quiet = TRUE)
  n <- length(metabolites_studied)
  for (number in 1:n) {
    while (!(metabolites_studied[number] %in% raw_data$metabolite)) {
      metabolite_not_found <- metabolites_studied[number]
      cat(paste("metabolite", metabolite_not_found, "not found, rewrite the list"))
      metabolites_studied <- scan(what = character(), quiet = TRUE)}
  }}



  #1.4 Select documents to include in graphs---------------------------
choice_number_document <- menu( 
  choices = c("Graphs for a single document", "Combined graphs for multiple documents"),
  title = "Please choose an option number:"
)
if (choice_number_document == 1) { #if one document select ask which file
  choice_file_selected <- menu(
    choices =  files,
    title = "Please choose an option number") 
} else {
  choice_file_selected <- 0}

if (choice_file_selected != 0){ #Filter by file if only one document was chosen
  selected_file <- files[choice_file_selected]
  raw_data <- raw_data |> 
    filter(file == selected_file)}

  #1.5 Restring injection dates to one per file-------------------------------------------
files <- unique(raw_data$file)

for (current_file in files) {
  file_rows <- raw_data$file == current_file
  injection_dates <- unique(raw_data$injection_date[file_rows])
  
  if (length(injection_dates) > 1) {
    cat(paste0("More than one injection date for file ", current_file, 
               ". Please select one:\n"))
    
    selected_index <- menu(injection_dates, title = "Select injection date:")
    selected_date <- injection_dates[selected_index]
    
    raw_data$injection_date[file_rows] <- selected_date
  }
}

#view(raw_metabolite)   
  #1.6 Create directories to save graphs --------------------------------

dir.create("Plasma/Comparisions/Young (12 weeks)/5mM 10mM 20mM IP-Glc", recursive = TRUE, showWarnings = FALSE)
dir.create("Plasma/Analytical control/Exact mass variation", recursive = TRUE, showWarnings = FALSE)
dir.create("Plasma/Analytical control/Retention time variation", recursive = TRUE, showWarnings = FALSE)



#ANALYTIC CONTROLS
#------------------------------------------------------------------------------  
#2. ANALYTIC CONTROLS-----------------------------------------------------
#2.1 EXACT MASS GRAPH--------------------------------------------------------------
for (metabolite_studied in metabolites_studied) {
 
  #2.1.1 Create database for each metabolite-----------------------------------------
   num_rows <- raw_data |> 
    filter(metabolite == metabolite_studied) |>
    nrow()
   
  raw_metabolite <- raw_data |>
    filter(metabolite == metabolite_studied)|>
    mutate("exact_mass" = runif(num_rows, min = 185, max = 195))|> #create random values
    select(sample, metabolite, exact_mass, isotopologue, injection_number, injection_date) 
  
  #2.1.2 Statistic tests------------------------------------------
  # Linear model
  model_lm <- lm(exact_mass ~ injection_number, data = raw_metabolite)
  lm_summary <- summary(model_lm)
  
  # GAM model
  model_gam <- gam(exact_mass ~ s(injection_number), data = raw_metabolite)
  gam_summary <- summary(model_gam)
  
  lm_pval <- coef(summary(model_lm))["injection_number", "Pr(>|t|)"]
  lm_rsq <- lm_summary$r.squared
  
  gam_pval <- gam_summary$s.table[1, "p-value"]  # p-value for s(injection_number)
  gam_edf <- gam_summary$s.table[1, "edf"]      # effective degrees of freedom
  gam_rsq <- gam_summary$r.sq
  
  #2.1.3 Annotation_text-------------------
  annotation_text <- paste0(
    "Linear model:\n",
    "p = ", signif(lm_pval, 3), 
    ", R² = ", signif(lm_rsq, 3), "\n",
    "GAM model:\n",
    "p = ", signif(gam_pval, 3),
    ", EDF = ", signif(gam_edf, 3), 
    ", R² = ", signif(gam_rsq, 3)
  )

  #2.1.4 Exact mass graph all isotopologues (em_graph_all)------------------------------------------
  em_graph_all <- ggplot(
    data=raw_metabolite,
    mapping = aes(x= injection_number, y = exact_mass)) +
    theme_classic() + 
    ylim(150, 515) +
    scale_x_continuous(breaks = unique(raw_metabolite$injection_number))+
    
    geom_point( aes(colour = isotopologue)) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = "yellow", size = 1) +
    geom_point()  +
    annotate("text", x = 17, 
             y = max(raw_metabolite$exact_mass) + 100, 
             label = annotation_text, hjust = 0, vjust = 1, size = 3.5) +
    labs(
      title= paste("Exact mass graph variation \nInjection date:", injection_date_chosen,  "\nMolecule studied:", metabolite_studied),
      subtitle = "Red = Linear model; Yellow = Nonlinear GAM smoother",
      x = "Injection number in the sequence",
      y = "Exact mass")
  
  em_graph_all #visualise graph
  
  #2.1.5 Exact mass graph each isotopologue separated (em_graph_separated)------------------------------------------
    em_graph_separated <- ggplot(
    raw_metabolite,
    aes(x = injection_number, y = exact_mass)) +
    theme_classic() + 
    ylim(150, 515) +
    scale_x_continuous(breaks = seq(0, max(raw_metabolite$injection_number), by = 5))+
    
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = "yellow", size = 1) +
    annotate("text", x = 3, 
             y = max(raw_metabolite$exact_mass) + 200, 
             label = annotation_text, hjust = 0, vjust = 1, size = 3.5) +
    geom_point( aes(colour = isotopologue)) +
    theme_classic() +
    labs(
      title = paste("Exact mass graph variation \nInjection date:", injection_date_chosen,  "\nMolecule studied:", metabolite_studied),
      subtitle = "Red = Linear model; Yellow = Nonlinear GAM smoother",
      x = "Injection number in the sequence",
      y = "Exact mass"
    ) + facet_wrap(
      vars(isotopologue),
      axes = "all",
      axis.labels = "all_y", 
      ncol=3,
      scales = "fixed"
    )
  
  em_graph_separated
  
  
  
  #2.1.6 Save exact mass graphs------------------------------------------------------------- 
  
  name_em_graph_a <- paste(injection_date_chosen, "_", metabolite_studied, "_exact_mass_all.png", sep = "")
  name_em_graph_s <- paste(injection_date_chosen, "_", metabolite_studied, "_exact_mass_separated.png", sep = "")
  ggsave(plot = em_graph_all,
         filename = file.path("Plasma", "Analytical control", "Exact mass variation", name_em_graph_a), 
         device = "png",
         width = 10,
         height = 7,
         units = "in",
         dpi = 300,
         limitsize = FALSE)
 
  ggsave(plot = em_graph_separated,
         filename = file.path("Plasma", "Analytical control", "Exact mass variation", name_em_graph_s), 
         device = "png",
         width = 10,
         height = 7,
         units = "in",
         dpi = 300,
         limitsize = FALSE)
  
}

#------------------------------------------------------------------------------
#2.2. RETENTION TIME GRAPH----------------------------------------------------------
for (metabolite_studied in metabolites_studied) {
  
  #2.2.1 Create database with reduced information--------------------------- 
  num_rows <- raw_data |> 
    filter(metabolite == metabolite_studied) |>
    nrow()
  
  raw_metabolite <- raw_data |>
    filter(metabolite == metabolite_studied)|>
    mutate("retention_time" = runif(num_rows, min = 6.54, max=7.54))|>
    select(sample, metabolite, isotopologue, retention_time, injection_number, injection_date) 
  #view(raw_metabolite) 
  
  
  #2.2.2 Statistic tests------------------------------------------
  # Linear model 
  model_lm <- lm(retention_time ~ injection_number, data = raw_metabolite)
  lm_summary <- summary(model_lm)
  
  # GAM model
  model_gam <- gam(retention_time ~ s(injection_number), data = raw_metabolite)
  gam_summary <- summary(model_gam)
  
  lm_pval <- coef(summary(model_lm))["injection_number", "Pr(>|t|)"]
  lm_rsq <- lm_summary$r.squared
  
  gam_pval <- gam_summary$s.table[1, "p-value"]  # p-value for s(injection_number)
  gam_edf <- gam_summary$s.table[1, "edf"]      # effective degrees of freedom
  gam_rsq <- gam_summary$r.sq
  
  #2.2.3 Annotation_text-------------------
  annotation_text <- paste0(
    "Linear model:\n",
    "p = ", signif(lm_pval, 3), 
    ", R² = ", signif(lm_rsq, 3), "\n",
    "GAM model:\n",
    "p = ", signif(gam_pval, 3),
    ", EDF = ", signif(gam_edf, 3), 
    ", R² = ", signif(gam_rsq, 3)
  )
  
  #2.2.4 Retention time graph all isotopologues (rt_graph_all)------------------------------------------
  rt_graph_all <- ggplot(
    data = raw_metabolite,
    mapping = aes(x = injection_number, y = retention_time)
  ) +
    theme_classic() +
    ylim(0, 12) +  
    scale_x_continuous(breaks = unique(raw_metabolite$injection_number)) +
    
    geom_point(aes(colour = isotopologue)) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = "yellow", size = 1) +
    annotate(
      "text", x = 17, 
      y = max(raw_metabolite$retention_time) + 2,  
      label = annotation_text, hjust = 0, vjust = 1, size = 3.5
    ) +
    labs(
      title = paste("Retention time graph variation \nInjection date:", injection_date_chosen, "\nMolecule studied:", metabolite_studied),
      subtitle = "Red = Linear model; Yellow = Nonlinear GAM smoother",
      x = "Injection number in the sequence",
      y = "Retention time"
    )
  
  rt_graph_all 
  
  #2.2.5 Retention time graph each isotopologue separated (rt_graph_separated)------------------------------------------
  rt_graph_separated <- ggplot(
    raw_metabolite,
    aes(x = injection_number, y = retention_time)
  ) +
    theme_classic() +
    ylim(0, 12) +  
    scale_x_continuous(breaks = seq(0, max(raw_metabolite$injection_number), by = 5))+
    
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = "yellow", size = 1) +
    annotate(
      "text", x = 3, 
      y = max(raw_metabolite$retention_time) + 4,  
      label = annotation_text, hjust = 0, vjust = 1, size = 3.5
    ) +
    geom_point(aes(colour = isotopologue)) +
    labs(
      title = paste("Retention time graph variation \nInjection date:", injection_date_chosen, "\nMolecule studied:", metabolite_studied),
      subtitle = "Red = Linear model; Yellow = Nonlinear GAM smoother",
      x = "Injection number in the sequence",
      y = "Retention time"
    ) +
    facet_wrap(
      vars(isotopologue),
      ncol = 3,
      scales = "fixed"
    )
  
  rt_graph_separated 
  
  #2.2.6 Save retention time graphs-------------------------------------------------------------
  
  name_rt_graph_a <- paste(injection_date_chosen, "_", metabolite_studied, "_retention_time_all.png", sep = "")
  name_rt_graph_s <- paste(injection_date_chosen, "_", metabolite_studied, "_retention_time_separated.png", sep = "")
  
  ggsave(
    plot = rt_graph_all,
    filename = file.path("Plasma", "Analytical control", "Retention time variation", name_rt_graph_a),
    device = "png",
    width = 10,
    height = 7,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  ggsave(
    plot = rt_graph_separated,
    filename = file.path("Plasma", "Analytical control", "Retention time variation", name_rt_graph_s),
    device = "png",
    width = 10,
    height = 7,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
}

#-------------------------------------------------------------------------------
#3. STATISTICAL COMPARISONS------------------------------------------------------

for (metabolite_studied in metabolites_studied){
  #3.1 Filter database for selected metabolite (metabolite_summary)--------------------------------------------------
  metabolite_summary <- raw_data|>
    filter(metabolite == metabolite_studied) |> #Keep only selected metabolite
    drop_na(glucose_condition)|> #Remove entries without infused or injected glucose
    arrange(glucose_condition, isotopologue) |> #sort by glucose concentration and isotopologue
    select(file,sample, metabolite, glucose_condition , isotopologue, isotopologue_fraction) |> 
    #only show reduced number of columns
    mutate("isotopologue_fraction_%" = isotopologue_fraction * 100)|> 
    select(-isotopologue_fraction)
  
  
  
  #view(metabolite_summary)
  
  #3.2 Keep only the most labeled isotopologue -----------------------
  highest_isotopologue <- max(metabolite_summary$isotopologue) #search for maximum
  metabolite_summary <- metabolite_summary |>
    filter(isotopologue == highest_isotopologue) 
  
  
  #3.3 Calculate mean and variance for each glucose condition------------------------
  metabolite_summary <- metabolite_summary |>
    group_by(glucose_condition )|>
    mutate(
      mean_i_f = mean(`isotopologue_fraction_%`),  
      sd_i_f=(sd(`isotopologue_fraction_%`)))
  
  #view(metabolite_summary)
  #3.4 Calculate biological replica (n) per condition--------------------------------------
  metabolite_summary1 <- metabolite_summary |>
    add_count(glucose_condition, name = "n")|>  
    group_by(glucose_condition)|>
    slice_head(n = 1)|> #Keep one row per condition
    select(glucose_condition, mean_i_f, sd_i_f, n,) 
  #view(metabolite_summary1) #show reduced version of the important information
  
  
  #3.5 Statistic test, create comparisons--------------------------------------------------------------
  
  expected_conditions <- c("5", "10", "20", "IP-Glc", "IP-Lac", "IP-Ala", "IP-BHB")
  available_conditions <- intersect(expected_conditions, unique(metabolite_summary$glucose_condition))
  if (!("IP" %in% available_conditions)) {
    comparisons_conditions <- combn(available_conditions, 2, simplify = FALSE)
  } else {
    group1 <- c("5", "10", "20")
    group2 <- intersect(expected_conditions, c("IP-Glc", "IP-Lac", "IP-Ala", "IP-BHB") ) 

  
  ip_comparisons <- unlist(lapply(group1, function(g1) {
    lapply(group2, function(g2) c(g1, g2))
  }), recursive = FALSE)
  
  glucose_comparisons <- combn(group1, 2, simplify = FALSE)
  
  comparisons_conditions <- c(glucose_comparisons, ip_comparisons)}
  
  
  
  #3.6 Create the plot(if_graph)---------------------------
  file_graph <- paste(unique(metabolite_summary$file), collapse = "\n")
  
  "if_graph" <-ggplot(metabolite_summary, aes(x = glucose_condition, y=`isotopologue_fraction_%`)) +
    theme_classic() +
    #scale_y_continuous(limits = c(0, 30)) + #fixed scale
    geom_bar(   #mean of isotopologue fraction for each glucose concentration
      data = metabolite_summary1,
      aes(x=glucose_condition, y=mean_i_f),
      stat = "identity", #if not R uses "count" as y axis
      alpha = 0.2,
      width = 0.4,
      inherit.aes = FALSE #if not aes and data used are the ones specified in ggplot
    ) + 
    
    geom_errorbar(      #add bars for standard deviation
      data = metabolite_summary1,
      aes(x = glucose_condition, ymin = mean_i_f - sd_i_f, ymax = mean_i_f + sd_i_f),
      width = 0.1,
      alpha=0.4,
      inherit.aes = FALSE) +
    
    geom_jitter(
      data = metabolite_summary,
      aes(x = glucose_condition, y = `isotopologue_fraction_%`, colour = file, shape = file),
      stroke = 1, #density of bords
      fill = "white",
      inherit.aes = FALSE,
      width = 0.05,
      height = 0)+
    
    theme(  #add legend
      legend.text = element_text(size = 10),
      legend.position = "bottom")+
    
    geom_text(                    #add mean and n next to the bar
      data = metabolite_summary1,
      aes(x =glucose_condition, y = mean_i_f +1,      #add mean and n in top of the sd bar
          label = paste("Mean: ", round(mean_i_f, 2),"% \n     n =", n)),
      inherit.aes = FALSE,
      size = 2.5,
      vjust= -0.1,        #adjust position x axis
      hjust= -0.3         #adjust position y axis
    ) +
    
    geom_signif(
      comparisons = comparisons_conditions,
      test = "t.test", #select statistic test
      map_signif_level = TRUE, #adds significance stars
      textsize = 2,
      margin_top = 0.03, #adds space to error bar
      step_increase = 0.1, #adds space between bars
      tip_length = 0.01) + #length of the vertical bit
    
    labs(  #add results of the tests
      title = paste("Isotopologue fractions in relation to \nconcentration of uniformly [U-¹³C] labeled glucose"),
      subtitle =paste("Metabolite studied:", metabolite_studied, "\n File shown:", file_graph),
      x = "Glucose concentration (mM) or injection",
      y= "Isotopologue fraction (%)") 
  
  
  if_graph #visualize plot
  
  
  #3.7 Save if_graph -----------------------------------------------------------
  injection_date_name <- paste(unique(raw_metabolite$injection_date), collapse = "_")
  name_if_graph <- paste(injection_date_name, "_", metabolite_studied, ".png", sep= "")
  ggsave(plot = if_graph,
         filename = file.path("Plasma", "Comparisions","Young (12 weeks)", "5mM 10mM 20mM IP-Glc", name_if_graph), 
         device = "png",
         width = 10,
         height = 6,
         units = "in",
         dpi = 300,
         limitsize = FALSE) 
} 


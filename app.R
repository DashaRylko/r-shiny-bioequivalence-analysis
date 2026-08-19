
## Полный и точный код приложения (app.R)

library(shiny)
library(bslib)
library(tidyverse)
library(readxl)
library(NonCompart)
library(ggplot2)
library(sasLM)
library(DT)
library(PowerTOST)
library(shinyWidgets)

options(shiny.maxRequestSize = 50*1024^2) # 50 MB
source('functions.R') 

# ==============================================================================
# ФУНКЦИИ ИЗ ОРИГИНАЛЬНЫХ ФАЙЛОВ (БЕЗ ИЗМЕНЕНИЙ В ЛОГИКЕ)
# ==============================================================================

nca_calculation <- function(data) {
  min_time <- data %>% select(Time) %>% min()
  data <- data %>%
    mutate(
      Concentration = ifelse(Concentration == 'BLOQ', NA, Concentration),
      across(c(Time, Concentration), as.numeric))
  
  data_subj_period <- data %>% 
    mutate(Subject = str_c(Subject,'_', Period), 
           Time = ifelse(Time==min_time,0,Time)
           )
  
  nca_results <- tblNCA(data_subj_period , key = "Subject", colTime = "Time", colConc =
                          "Concentration", down = "Log") %>% 
    select(Subject, CMAX, TMAX, CLST, TLST, CLSTP, LAMZNPT, LAMZHL, LAMZ, AUCLST, AUCIFP)
  
  nca_results_with_form <- nca_results %>%
    left_join(data_subj_period %>% select(Subject, Formulation,Sequence) %>% unique(),
              by = 'Subject') %>%
    separate(Subject, into = c('Subject', 'Period'), sep = '_')
  return(nca_results_with_form)
}

nca_rounding <- function(nca_table, time_rounding, pk_rounding) {
  nca_table %>% 
    mutate(across(c(TMAX,TLST,LAMZHL),~round(.,time_rounding)),
           across(c(CMAX,CLST,CLSTP,AUCLST,AUCIFP),~round(.,pk_rounding)),
           LAMZ = LAMZ %>% round(4),
           LAMZNPT = LAMZNPT %>% round(0))
}

individ_rounding <- function(indidvid_table, rounding) {
  indidvid_table %>% 
    mutate(across(where(is.numeric),~round(.,rounding)))
}

data_individ_conc_calc <- function(data) {
  data_individ_conc <- data %>%
    mutate(Concentration = Concentration %>% as.numeric()) %>%
    pivot_wider(names_from = Time, values_from = Concentration) %>%
    mutate(Subject = paste0(Subject, ' П', Period, ' ', Formulation)) %>%
    select(-c(Sequence, Formulation))
  return(data_individ_conc)
}

desc_statistic <- function(table) {
  statistics <- list(
    `__N` = ~sum(!is.na(.x)),
    `__Mean` = ~ifelse(sum(!is.na(.x)) == 0, NA, mean(.x, na.rm = TRUE)),
    `__SD` = ~ifelse(sum(!is.na(.x)) == 0, NA,sd(.x, na.rm = TRUE)),
    `__Median` = ~ifelse(sum(!is.na(.x)) == 0, "-", median(.x, na.rm = TRUE)),
    `__Min` = ~ifelse(sum(!is.na(.x)) == 0, "-", min(.x, na.rm = TRUE)),
    `__Max` = ~ifelse(sum(!is.na(.x)) == 0, "-", max(.x, na.rm = TRUE)),
    `__CV%` = ~ifelse(sum(!is.na(.x)) == 0, "-", (sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)) * 100),
    `__GeomMean` = ~ifelse(sum(!is.na(.x)) == 0,"-", exp(sum(log(.x[.x > 0]), na.rm = TRUE) / length(.x[.x > 0])))
  )
  table %>%
    summarise(across(where(is.numeric) & !c(Period), statistics)) %>% 
    pivot_longer(cols = everything(),names_to = 'name', values_to = 'value') %>% 
    separate(name, into = c("Variable", "Statistics"), sep="___") %>% 
    pivot_wider(names_from = "Variable", values_from = "value")
}

# --- Логика BE из cross.txt ---
be_result_cross <- function(PK = 'CMAX', nca_data) {
  f1 = log(nca_data[[PK]]) ~ Sequence / Subject + Period + Formulation
  be_res <- list(
    anova = as.data.frame(unclass(GLM(f1, nca_data)$`ANOVA`)) %>% mutate(across(where(is.numeric), ~ round(.,5))),
    type3 = as.data.frame(unclass(GLM(f1, nca_data)$`Type III`)) %>% mutate(across(where(is.numeric), ~ round(.,5))),
    seqtest = as.data.frame(unclass(T3test(f1, nca_data, H = "Sequence", E = c("Sequence:Subject"))$Sequence)) %>% mutate(across(where(is.numeric), ~ round(.,5))),
    CIres = cbind(Measure = c('normal, %','ln'), as.data.frame(unclass(rbind(exp(CIest(f1, nca_data, 'Formulation', c(-1, 1), 0.9)), CIest(f1, nca_data, 'Formulation', c(-1, 1), 0.9))))) %>%
      select('Measure','Estimate', 'Lower CL', 'Upper CL') %>%
      mutate(across(where(is.numeric), ~ if_else(Measure == 'normal, %', round(. * 100, 2) %>% as.character(), round(., 5) %>% as.character())))
  )
  be_res$anova_total = rbind(be_res$seqtest) %>% rbind(be_res$type3) %>% rbind(be_res$anova) %>%
    rownames_to_column() %>% 
    filter(rowname %in% c('Sequence', 'Sequence:Subject', 'Period', 'Formulation', 'RESIDUALS')) %>% 
    rename('Parameter'= rowname)
  be_res$cv = 100*sqrt(exp(be_res$anova['RESIDUALS','Mean Sq'])-1)
  be_res$power = power.TOST(
    CV = be_res$cv/100, 
    theta0 = exp(CIest(f1, nca_data, 'Formulation', c(-1, 1), 0.9))[,'Estimate'], 
    n = nca_data %>% group_by(Sequence) %>% summarise(n = n()) %>% pull(n)/2,
    design = "2x2", method="exact"
  )
  return(be_res)
}

# --- Логика BE из replicate.txt ---
be_result_replicate <- function(PK = 'CMAX', nca_data, alpha = 0.05, design = '2x2x4') {
  nca_data <- nca_data %>% 
    mutate(Subject = parse_number(Subject),
           Sequence = case_when(Sequence == 1 ~ 'TRTR', Sequence == 2 ~ 'RTRT', TRUE ~ as.character(Sequence)),
           Period = Period %>% as.numeric(),
           Sequence = Sequence %>% as.factor()) %>% 
    rename('PK' := {{PK}},'subject'=Subject,'period'=Period,'treatment'=Formulation,'sequence'=Sequence) %>% 
    select(PK,subject,period, treatment, sequence)
  
  nca_data$subject <- as.factor(as.character(nca_data$subject))
  nca_data$period <- as.factor(as.character(nca_data$period))
  nca_data$sequence <- as.factor(as.character(nca_data$sequence))
  nca_data$treatment <- as.factor(as.character(nca_data$treatment))
  
  res <- method.A(data = nca_data, print = FALSE, ext = '', alpha = alpha, verbose = TRUE, file = "project_name", path.in = "")
  
  be_res <- list(
    anova = as.data.frame(res$anova) %>% mutate(across(where(is.numeric), ~ round(.,5))),
    CIres = cbind(Measure = c('normal, %','ln', 'BE borders, %'), as.data.frame(unclass(rbind(
      c(res$res$`PE(%)`,res$res$`CL.lo(%)`,res$res$`CL.hi(%)`),
      c(log(res$res$`PE(%)`/100),log(res$res$`CL.lo(%)`/100),log(res$res$`CL.hi(%)`/100)),
      c(res$res$`PE(%)`,res$res$`L(%)`,res$res$`U(%)`)
    )))) %>%
      rename('Estimate' = 'V1', 'Lower CL' = 'V2', 'Upper CL' = 'V3') %>%
      mutate(across(where(is.numeric), ~ if_else(Measure == 'normal, %' | Measure == 'BE borders, %', round(., 2) %>% as.character(), round(., 5) %>% as.character())))
  )
  be_res$anova = be_res$anova %>% rownames_to_column() %>% 
    filter(rowname %in% c('sequence', 'sequence:subject', 'period', 'treatment', 'Residuals')) %>% 
    rename('Parameter'= rowname)
  be_res$cv_r = res$res$`CVwR(%)`
  be_res$cv_t = res$res$`CVwT(%)`
  be_res$power = power.TOST(alpha = alpha, CV = be_res$cv_r/100, theta0 = res$res$`PE(%)`/100, n = nca_data %>% group_by(sequence) %>% summarise(n = n()) %>% pull(n)/4, design = design, method="exact")
  return(be_res)
}

# --- Функции графиков из cross.txt ---
single_subject_plot_cross <- function(subj_data, y_coord = c('norm','ln'), max_time, breaks) {
  subj_id <- unique(subj_data$Subject)
  norm_plot <- subj_data %>% ggplot(aes(x = Time, y = Concentration, color = Formulation, group = Formulation)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4")) +
    scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = paste0(subj_id), x = "Время (ч)", y = "Концентрация (нг/мл)", color = "Препарат:")
  
  ln_plot <- subj_data %>% ggplot(aes(x = Time, y = log(Concentration), color = Formulation, group = Formulation)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4")) +
    scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = paste0(subj_id), x = "Время (ч)", y = "Логарифм концентрации (нг/мл)", color = "Препарат:")
  if(y_coord == 'norm') { return(norm_plot) } else {return(ln_plot)}
}

all_subject_plots_cross <- function(data, y_coord = c('norm','ln'),max_time, breaks) {
  min_time <- data %>% select(Time) %>% min()
  data <- data %>% mutate(Concentration = ifelse(Concentration == 'BLOQ', NA, Concentration), across(c(Time, Concentration), as.numeric), Time = ifelse(Time==min_time,0,Time))
  plot_list <- lapply(unique(data$Subject), function(s) single_subject_plot_cross(data %>% filter(Subject == s), y_coord = y_coord, max_time, breaks))
  return(plot_list)
}

mean_subject_plots_cross <- function(data, y_coord = c('norm','ln'), max_time, breaks) {
  data <- data %>% mutate(Concentration = ifelse(Concentration == 'BLOQ', NA, Concentration), across(c(Time, Concentration), as.numeric))
  min_time <- data %>% select(Time) %>% min()
  data <- data %>% mutate(Time = ifelse(Time==min_time,0,Time)) %>% group_by(Time, Formulation) %>% summarise(Concentration = mean(Concentration, na.rm = T),.groups='drop')
  
  norm_plot <- data %>% ggplot(aes(x = Time, y = Concentration, color = Formulation, group = Formulation)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4")) +
    scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Усредненный профиль', x = "Время (ч)", y = "Концентрация (нг/мл)", color = "Препарат:")
  
  ln_plot <- data %>% ggplot(aes(x = Time, y = log(Concentration), color = Formulation, group = Formulation)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4")) +
    scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Усредненный профиль', x = "Время (ч)", y = "Логарифм концентрации (нг/мл)", color = "Препарат:")
  
  if(y_coord == 'norm') { return(norm_plot) } else {return(ln_plot)}
}
overlap_plot_cross <- function(data, y_coord = c('norm','ln'), form = c('T','R'), max_time, breaks) {
  data <- data %>% mutate(Concentration = ifelse(Concentration == 'BLOQ', NA, Concentration), across(c(Time, Concentration), as.numeric))
  min_time <- data %>% select(Time) %>% min()
  data <- data %>% mutate(Time = ifelse(Time==min_time,0,Time))
  norm_plot <- data %>% filter(Formulation == form) %>% ggplot(aes(x = Time, y = Concentration, group = Subject)) +
    geom_line(linewidth = 0.9, alpha = 1, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4")) +
    scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Усредненный профиль', x = "Время (ч)", y = "Концентрация (нг/мл)")
  ln_plot <- data %>% filter(Formulation == form) %>% ggplot(aes(x = Time, y = log(Concentration), group = Subject)) +
    geom_line(linewidth = 0.9, alpha = 1, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4")) +
    scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Усредненный профиль', x = "Время (ч)", y = "Логарифм концентрация (нг/мл)")
  if(y_coord == 'norm') { return(norm_plot) } else {return(ln_plot)}
}
## --- Функции графиков из replicate.txt ---
single_subject_plot_rep <- function(subj_data, y_coord = c('norm','ln'), max_time, breaks) {
  subj_id <- unique(subj_data$Subject)
  norm_plot <- subj_data %>% ggplot(aes(x = Time, y = Concentration, color = Legend_Group, group = Legend_Group)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = paste0(subj_id), x = "Время (ч)", y = "Концентрация (нг/мл)", color = '')
  ln_plot <- subj_data %>% ggplot(aes(x = Time, y = log(Concentration), color = Legend_Group, group = Legend_Group)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_x_continuous(breaks = seq(0,max_time, max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = paste0(subj_id), x = "Время (ч)", y = "Логарифм концентрации (нг/мл)", color = '')
  if(y_coord == 'norm') { return(norm_plot) } else {return(ln_plot)}
}
all_subject_plots_rep <- function(data, y_coord = c('norm','ln'), max_time, breaks) {
  min_time <- data %>% select(Time) %>% min()
  data <- data %>% mutate(Legend_Group = paste0("Период: ", Period, " Препарат: ", Formulation), Period = Period %>% as.factor(), Concentration = ifelse(Concentration == 'BLOQ', NA, Concentration), across(c(Time, Concentration), as.numeric), Time = ifelse(Time==min_time,0,Time))
  plot_list <- lapply(unique(data$Subject), function(s) single_subject_plot_rep(data %>% filter(Subject == s), y_coord = y_coord, max_time, breaks))
  return(plot_list)
}
mean_subject_plots_rep <- function(data, y_coord = c('norm','ln'), max_time, breaks) {
  data <- data %>% mutate(Concentration = ifelse(Concentration == 'BLOQ', NA, Concentration), across(c(Time, Concentration), as.numeric))
  min_time <- data %>% select(Time) %>% min()
  data <- data %>% mutate(Time = ifelse(Time==min_time,0,Time)) %>% group_by(Time, Formulation) %>% summarise(Concentration = mean(Concentration, na.rm = T),.groups='drop')
  norm_plot <- data %>% ggplot(aes(x = Time, y = Concentration, color = Formulation, group = Formulation)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4"))+
    scale_x_continuous(breaks = seq(0,max_time,max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Усредненный профиль', x = "Время (ч)", y = "Концентрация (нг/мл)", color = "Препарат:")
  ln_plot <- data %>% ggplot(aes(x = Time, y = log(Concentration), color = Formulation, group = Formulation)) +
    geom_line(linewidth = 1.05, alpha = 0.8, na.rm = T) + scale_color_manual(values = c("T" = "#c7000d", "R" = "#2278d4"))+
    scale_x_continuous(breaks = seq(0,max_time,max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Усредненный профиль', x = "Время (ч)", y = "Логарифм концентрации (нг/мл)", color = "Препарат:")
  if(y_coord == 'norm') { return(norm_plot) } else {return(ln_plot)}
}
overlap_plot_rep <- function(data, y_coord = c('norm','ln'), form = c('T','R'), max_time, breaks) {
  data <- data %>% mutate(Concentration = ifelse(Concentration == 'BLOQ', NA, Concentration), across(c(Time, Concentration), as.numeric))
  min_time <- data %>% select(Time) %>% min()
  data <- data %>% mutate(Time = ifelse(Time==min_time,0,Time))
  norm_plot <- data %>% filter(Formulation == form) %>% ggplot(aes(x = Time, y = Concentration, group = interaction(Subject,Period))) +
    geom_line(linewidth = 0.6, alpha = 1, na.rm = T) + scale_x_continuous(breaks = seq(0,max_time,max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Профиль в наложении', x = "Время (ч)", y = "Концентрация (нг/мл)")
  ln_plot <- data %>% filter(Formulation == form) %>% ggplot(aes(x = Time, y = log(Concentration), group = interaction(Subject,Period))) +
    geom_line(linewidth = 0.6, alpha = 1, na.rm = T) + scale_x_continuous(breaks = seq(0,max_time,max_time/breaks), limits = c(0,max_time)) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) +
    labs(title = 'Профиль в наложении', x = "Время (ч)", y = "Логарифм концентрация (нг/мл)")
  if(y_coord == 'norm') { return(norm_plot) } else {return(ln_plot)}
}
## ==============================================================================## 2. ИНТЕРФЕЙС ПРИЛОЖЕНИЯ (UI) С ОРИГИНАЛЬНЫМ СТИЛЕМ СSS## ==============================================================================
# Custom CSS - единый стиль для всех окон
custom_css <- "
/* Навигационная панель */
.navbar {
  background-color: #6fb856 !important;
  padding: 0 20px;
  height: 60px;
  display: flex;
  align-items: center;
}

/* Заголовок */
.navbar-brand {
  color: white !important;
  font-weight: bold;
  font-size: 1.5rem;
  padding: 0;
  margin-right: 40px;
  height: 60px;
  display: flex;
  align-items: center;
}

/* Контейнер вкладок */
.navbar-nav {
  display: flex;
  align-items: center;
  height: 60px;
  margin: 0;
  padding: 0;
}

/* Вкладки */
.nav-item {
  height: 60px;
  display: flex;
  align-items: center;
}

.nav-link {
  color: rgba(255, 255, 255, 0.85) !important;
  font-weight: 500;
  font-size: 1rem;
  padding: 0 20px !important;
  height: 60px;
  display: flex;
  align-items: center;
  transition: all 0.2s;
  border-bottom: 3px solid transparent;
  margin: 0 2px;
}

.nav-link:hover {
  color: white !important;
  background-color: rgba(255, 255, 255, 0.1) !important;
  border-bottom: 3px solid rgba(255, 255, 255, 0.3);
}

.nav-link.active {
  color: white !important;
  font-weight: 600;
  background-color: rgba(255, 255, 255, 0.15) !important;
  border-bottom: 3px solid white;
}

/* Основной контент */
.main-content {
  padding-top: 20px;
}

/* Стиль карточек как в Example Files */
.unified-card {
  border: 1px solid #d1e3ff !important;
  border-radius: 8px !important;
  background-color: #f8fbfe !important;
  margin-bottom: 20px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.05);
}

.unified-card .card-header {
  background-color: #e8f1ff !important;
  border-bottom: 1px solid #d1e3ff !important;
  font-weight: 600;
  color: #6fb856;
  border-radius: 8px 8px 0 0 !important;
  padding: 15px 20px;
}

.unified-card .card-body {
  padding: 20px;
}

/* Карточки с цветными акцентами */
.card-accent-blue {
  border-left: 4px solid #aedba0 !important;
}

.card-accent-green {
  border-left: 4px solid #aedba0 !important;
}

.card-accent-red {
  border-left: 4px solid #aedba0 !important;
}

.card-accent-orange {
  border-left: 4px solid #aedba0 !important;
}

.card-accent-purple {
  border-left: 4px solid #aedba0 !important;
}

/* Кнопки в едином стиле */
.unified-btn {
  background-color: #6fb856;
  color: white;
  border: none;
  padding: 10px 20px;
  border-radius: 4px;
  font-weight: 500;
  transition: all 0.2s;
  width: 100%;
  margin-top: 10px;
}

.unified-btn:hover {
  background-color: #0eb00e;
  color: white;
  transform: translateY(-1px);
  box-shadow: 0 4px 8px rgba(0,0,0,0.1);
}

.unified-btn-primary {
  background-color: #6fb856;
}

.unified-btn-success {
  background-color: #6fb856;
}

.unified-btn-success:hover {
  background-color: #0eb00e;
}

/* Поля ввода */
.unified-input {
  border: 1px solid #d1e3ff;
  border-radius: 4px;
  padding: 8px 12px;
  font-size: 0.95rem;
}

.unified-input:focus {
  border-color: #6fb856;
  box-shadow: 0 0 0 0.2rem rgba(30, 75, 140, 0.25);
}

/* Таблицы */
.unified-table {
  border: 1px solid #d1e3ff;
  border-radius: 8px;
  overflow: hidden;
}

.unified-table .dataTables_wrapper {
  border-radius: 8px;
}

/* Секция с примерами */
.example-section {
  background-color: #f8fbfe;
  border: 1px solid #d1e3ff;
  border-radius: 8px;
  margin-bottom: 20px;
}

.example-section .card-header {
  background-color: #e8f1ff;
  border-bottom: 1px solid #d1e3ff;
}

/* Уведомления в едином стиле */
.unified-alert {
  border-radius: 6px;
  padding: 12px 15px;
  margin: 10px 0;
  border: 1px solid transparent;
}

.unified-alert-success {
  background-color: #d4edda;
  border-color: #c3e6cb;
  color: #155724;
}

.unified-alert-danger {
  background-color: #f8d7da;
  border-color: #f5c6cb;
  color: #721c24;
}

.unified-alert-info {
  background-color: #d1ecf1;
  border-color: #bee5eb;
  color: #0c5460;
}

/* Иконки справа в навбаре */
.navbar-right {
  margin-left: auto;
  display: flex;
  align-items: center;
  height: 60px;
}

.navbar-right a {
  color: white !important;
  text-decoration: none;
  padding: 0 15px;
  display: flex;
  align-items: center;
}

.navbar-right a:hover {
  opacity: 0.8;
}

/* Поля ввода файлов */
.file-input-wrapper .form-control {
  height: 38px;
  border: 1px solid #d1e3ff;
  border-radius: 4px;
}

.file-input-wrapper .input-group-btn .btn {
  height: 38px;
  padding: 8px 15px;
  background-color: #6fb856;
  color: white;
  border: 1px solid #6fb856;
}

/* Слайдер */
.unified-slider .irs-bar {
  background-color: #6fb856;
}

.unified-slider .irs-handle {
  border: 3px solid #6fb856;
}

/* Селекторы */
.selectize-control.single .selectize-input {
  border: 1px solid #d1e3ff;
  border-radius: 4px;
}

.selectize-control.single .selectize-input:focus {
  border-color: #6fb856;
}

.btn-group .btn {
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
}
      
/* Вкладки с графиками */
/* Цвет неактивных вкладок */
.nav-tabs .nav-link { 
  color: #6fb856 !important; 
  padding: 0 15px;
  display: flex;
  align-items: center;
}       

/* Цвет активной вкладки */
.nav-tabs .nav-link.active { 
  color: #e74c3c !important; 
  padding: 0 15px;
  display: flex;
  align-items: center;
} 
"
ui <- fluidPage(
  tags$head(tags$style(HTML(custom_css))),
  navbarPage(
    title = div(style = "display: flex; align-items: center; height: 60px;", "BE analysis"),
    id = "main_navbar",
    theme = bs_theme(version = 5, primary = "#6fb856", secondary = "#ba5536", success = "#52b788", info = "#693d3d", warning = "#ffb703", danger = "#e63946", font_scale = 0.95),
    tabPanel(
      title = "Upload Your Data", icon = icon("cloud-upload-alt"),
      div(
        class = "main-content",
        # Секция с примерами файлов
        card(
          class = "unified-card example-section",
          card_header("Example Files - Download and Review Format"),
          div(
            style = "padding: 20px;",
            layout_columns(
              col_widths = c(4, 4, 4),
              height = "auto",
              # Карточка 1
              div(
                style = "height: 100%;",
                card(class = "unified-card card-accent-blue",
                     card_header("Choose Design of Study"),
                     card_body(
                       selectInput("design_type", label = "Choose the Design of Study:", choices = c("Cross-over 2x2" = "cross", "Replicative 2x2x4, 2x2x3, 2x3x3" = "replicate"), selected = "cross")
                     )
                ), 
                card(
                  class = "unified-card card-accent-blue",
                  card_header("Data File Examples"),
                  card_body(
                    p("Main dataset for ", strong("crossover design")," with concentration values.", class = "text-muted"),
                    downloadButton("download_example_data_cross", label = "Download CSV", class = "unified-btn"),
                    p("Main dataset for ",strong("replicate design")," with concentration values.", class = "text-muted"),
                    downloadButton("download_example_data_rep", label = "Download CSV", class = "unified-btn")
                  )
                )
              ),
              
              # Карточка 3
              div(
                style = "height: 100%;",
                card(class = "unified-card card-accent-red",
                     card_header("Format Requirements"),
                     style = "height: 100%;",
                     card_body(
                       layout_column_wrap(
                         tags$ul(
                           style = "padding-left: 20px; margin: 0;",
                           tags$li("xlsx/csv format with UTF-8 encoding"),
                           tags$li("No missing values in key fields"),
                           tags$li("Headers exactly as in examples")
                         ),
                         tags$ul(
                           style = "padding-left: 20px; margin: 0;",
                           tags$li("Subject format: 'Sxx'"),
                           tags$li("Time: decimal number (1.333 e.g.)"),
                           tags$li("Concentration: decimal number (20.536 e.g.)"),
                           tags$li("Formulation: 'R' or 'T' only"),
                           tags$li("Sequence: 'RT' or 'TR', 'TRTR' or 'RTRT' e.g."),
                           tags$li("Period: even number (1,2 e.g.)")
                           
                         )
                       )
                     ))
              ),
              
              # Карточка инструкций
              card(
                class = "unified-card card-accent-purple",
                card_header("Instructions"),
                card_body(
                  tags$ol(
                    style = "padding-left: 20px; margin: 0;",
                    tags$li("Download example files"),
                    tags$li("Prepare your data in the same format"),
                    tags$li("Upload file"),
                    tags$li("Click 'Load and Process Data'"),
                    tags$li("Wait for confirmation"),
                    tags$li("Go to 'NCA Results' or 'BE Results' tab")
                  ),
                  hr(style = "margin: 20px 0; border-color: #d1e3ff;"),
                  div(
                    style = "font-size: 0.9em; color: #666;",
                    icon("info-circle", style = "color: #6fb856;"),
                    " Maximum file size: 50 MB per file",
                    br(),
                    icon("info-circle", style = "color: #6fb856;"),
                    " Supported format: CSV and xlsx (UTF-8 encoding)"
                  )
                )
              )
              
              
            )
          )
        ),
        
        # Секция загрузки файлов и инструкций
        layout_columns(
          col_widths = c(6, 6),
          height = "auto",
          # Карточка загрузки файлов
          card(
            class = "unified-card card-accent-orange",
            card_header("Upload Your Files"),
            card_body(
              div(
                class = "file-input-wrapper",
                style = "margin-bottom: 20px;",
                fileInput(
                  "user_data",
                  label = div(icon("file-medical"), "Choose Main Data File (CSV, xlsx)"),
                  accept = c(".csv", "text/csv",".xlsx", ".xls"),
                  width = "100%",
                  buttonLabel = "Browse",
                  placeholder = "No file selected"
                )
              ),
              
              actionButton(
                "load_user_data",
                label = div(icon("play-circle"), "Load and Process Data"),
                class = "unified-btn unified-btn-success"
              ),
              
              uiOutput("upload_status"),
              
              div(
                class = "file-input-wrapper",
                style = "margin-bottom: 20px;",
                fileInput(
                  "user_data_mean",
                  label = div(icon("file-medical"), "Choose Data File For Mean Profile (CSV, xlsx)"),
                  accept = c(".csv", "text/csv",".xlsx", ".xls"),
                  width = "100%",
                  buttonLabel = "Browse",
                  placeholder = "No file selected"
                )
              ),
              
              actionButton(
                "load_user_data_mean",
                label = div(icon("play-circle"), "Load and Process Data"),
                class = "unified-btn unified-btn-success"
              ),
              
              uiOutput("upload_status_mean")
            )
          ),
          
          # Предпросмотр данных
          card(
            class = "unified-card card-accent-blue",
            card_header(
              div(
                style = "display: flex; justify-content: space-between; align-items: center;",
                span("Main Data Preview (First 10 rows)"),
                tags$small(textOutput("data_rows"), class = "text-muted")
              )
            ),
            div(class = "unified-table",
                style = "padding: 10px;",
                DTOutput("data_preview"))
          )
          
        )
      )
    ),
    tabPanel(
      title = 'NCA Results',
      icon = icon("table"),
      conditionalPanel(
        condition = "output.user_data_loaded == false",
        div(
          class = "main-content",
          style = "text-align: center; padding: 100px 20px;",
          card(
            class = "unified-card",
            style = "max-width: 600px; margin: 0 auto;",
            card_header("No Data Loaded"),
            card_body(
              icon("database", style = "font-size: 60px; color: #6c757d; margin-bottom: 20px;"),
              h4("Please upload your data first", style = "color: #6c757d;"),
              p("Go to the 'Upload Your Data' tab to load your BE data."),
              actionButton("go_to_upload",
                           label = div(icon("cloud-upload-alt"), "Go to Upload Data"),
                           class = "unified-btn",
                           style = "margin-top: 20px;"
              )
            )
          )
        )
      ),
      conditionalPanel(
        condition = "output.user_data_loaded == true",
        fluidRow(
          selectInput(
            inputId = 'tr',
            label = strong('Select formulation'),
            choices = c('T', 'R')
          ),
          numericInput(
            inputId = 'timeround',
            label = strong('Select time parameters rounding'),
            value = 3,
            min = 0,
            max = 10
          ),
          numericInput(
            inputId = 'pkround',
            label = strong('Select PK parameters rounding'),
            value = 3,
            min = 0,
            max = 10
          ),
        ),
        
        
        
        card(
          class = "unified-card card-accent-blue",
          card_header(
            div(
              style = "display: flex; justify-content: space-between; align-items: center;",
              span("NCA results")
            )
          ),
          
          div(class = "unified-table",
              style = "padding: 10px;",
              DTOutput('resultNCA')),
          
          div(class = "unified-table",
              style = "padding: 10px;",
              DTOutput('nca_desc_statistic'))
        )
      )
      
    ),
    tabPanel(
      title = 'BE Results',
      icon = icon("square-poll-vertical"),
      conditionalPanel(
        condition = "output.user_data_loaded == false",
        div(
          class = "main-content",
          style = "text-align: center; padding: 100px 20px;",
          card(
            class = "unified-card",
            style = "max-width: 600px; margin: 0 auto;",
            card_header("No Data Loaded"),
            card_body(
              icon("database", style = "font-size: 60px; color: #6c757d; margin-bottom: 20px;"),
              h4("Please upload your data first", style = "color: #6c757d;"),
              p("Go to the 'Upload Your Data' tab to load your BE data."),
              actionButton("go_to_upload",
                           label = div(icon("cloud-upload-alt"), "Go to Upload Data"),
                           class = "unified-btn",
                           style = "margin-top: 20px;"
              )
            )
          )
        )
      ),
      conditionalPanel(
        condition = "output.user_data_loaded == true",
        div(style = "padding: 20px;", uiOutput("dynamic_be_inputs")),
        card(
          class = "unified-card card-accent-blue",
          style = "margin: 0px 20px 20px 20px;",
          card_header(span("90 % confidence interval")),
          div(class = "unified-table", style = "padding: 10px;", DTOutput('resultCI'))
        ),
        card(
          class = "unified-card card-accent-blue",
          style = "margin: 0px 20px 20px 20px;",
          card_header(span("Anova total results")),
          div(class = "unified-table", style = "padding: 10px;", DTOutput('anova_total'))
        ),
        card(
          class = "unified-card card-accent-blue",
          style = "margin: 0px 20px 20px 20px;",
          card_header(span("CV and Power estmation")),
          div(class = "unified-table", style = "padding: 10px;", uiOutput('cv_power'))
        )
      )
    ), 
    tabPanel(
      title = 'Plots',
      icon = icon("chart-area"),
      conditionalPanel(
        condition = "output.user_data_loaded == false",
        div(
          class = "main-content",
          style = "text-align: center; padding: 100px 20px;",
          card(
            class = "unified-card",
            style = "max-width: 600px; margin: 0 auto;",
            card_header("No Data Loaded"),
            card_body(
              icon("database", style = "font-size: 60px; color: #6c757d; margin-bottom: 20px;"),
              h4("Please upload your data first", style = "color: #6c757d;"),
              p("Go to the 'Upload Your Data' tab to load your BE data."),
              actionButton("go_to_upload",
                           label = div(icon("cloud-upload-alt"), "Go to Upload Data"),
                           class = "unified-btn",
                           style = "margin-top: 20px;"
              )
            )
          )
        )
      ),
      conditionalPanel(condition = "output.user_data_loaded == true",
                       layout_sidebar(
                         sidebar = sidebar(
                           width = 320,
                           open = "always",
                           class = "unified-card",
                           style = "height: calc(100vh - 80px); overflow-y: auto;",
                           div(
                             class = "card-header",
                             h4("Data loading for mean profile", style = "margin: 0; color: #1e4b8c;")
                           ),

                           numericInput(inputId = 'max_time', label = 'Input max time', value = 72, min = 0, max = 1000000),
                           numericInput(inputId = 'breaks', label = 'Input number of breaks', value = 6, min = 0, max = 100)
                         ),
                         div(
                           class = "main-content",
                           navset_card_tab(
                             full_screen = TRUE,
                             title = "Plots",
                             nav_panel(
                               "Individual profiles",
                               card_title("Individual profiles"),
                               style = "padding: 20px;",
                               card(
                                 class = "unified-card card-accent-blue",
                                 card_header(span("Plots")),
                                 div(
                                   style = "display: flex; align-items: flex-end; gap: 15px; width: 100%;",
                                   
                                   # Внедряем кусочек CSS, который уберет скрытый отступ Bootstrap именно внутри этой группы
                                   tags$style(HTML("
    .no-margin-group .form-group { margin-bottom: 0 !important; }
    .no-margin-group .shiny-input-container { margin-bottom: 0 !important; }
  ")),
                                   
                                   div(
                                     class = "no-margin-group",
                                     style = "width: 30%;",
                                     radioGroupButtons(
                                       inputId = "y_coord",
                                       label = "Select option:",
                                       choices = c("Normal" = "norm", "Log" = "ln"),
                                       status = "primary",
                                       justified = TRUE,
                                       width = "100%"
                                     )
                                   ),
                                   
                                   div(
                                     style = "width: 30%;",
                                     downloadButton(
                                       outputId = "download_all_plots",
                                       label = "Скачать все индивидуальные графики (ZIP)",
                                       class = "btn-primary",
                                       style = "width: 100%; margin-bottom: 0 !important;"
                                     )
                                   )
                                 ),
                                 
                                 
                                 
                                 div(
                                   class = "unified-table",
                                   style = "padding: 10px;",
                                   uiOutput("plots_container")
                                 )
                               )
                             ),
                             nav_panel(
                               "Overlap profiles",
                               card_title("Overlap profiles"),
                               style = "padding: 20px;",
                               card(
                                 class = "unified-card card-accent-blue",
                                 card_header(span("Plots")),
                                 div(
                                   style = "display: flex; align-items: flex-end; gap: 15px; width: 100%;",
                                   
                                   # Внедряем кусочек CSS, который уберет скрытый отступ Bootstrap именно внутри этой группы
                                   tags$style(HTML("
    .no-margin-group .form-group { margin-bottom: 0 !important; }
    .no-margin-group .shiny-input-container { margin-bottom: 0 !important; }
  ")),
                                   
                                   div(
                                     class = "no-margin-group",
                                     style = "width: 30%;",
                                     radioGroupButtons(
                                       "y_coord_overlap",
                                       "Select option:",
                                       choices = c("Normal" = "norm", "Log" = "ln"),
                                       status = "primary",
                                       justified = TRUE,
                                       width = "100%"
                                     )
                                   ),
                                   
                                   div(
                                     class = "no-margin-group",
                                     style = "width: 30%;",
                                     radioGroupButtons(
                                       "formulation",
                                       "Select option:",
                                       choices = c("R" = "R", "T" = "T"),
                                       status = "primary",
                                       justified = TRUE,
                                       width = "100%"
                                     )
                                   ),
                                   
                                   div(
                                     style = "width: 30%;",
                                     downloadButton(
                                       outputId = "download_overlap_plot", 
                                       label = "Скачать профиль в наложении", 
                                       class = "btn-primary",
                                       style = "width: 100%; margin-bottom: 0 !important;"
                                     )
                                   )
                                 ),
                                 
                                 div(
                                   class = "unified-table",
                                   style = "padding: 10px;",
                                   plotOutput("overlap_plot")
                                 )
                               )
                             ),
                             nav_panel(
                               "Mean profile",
                               card_title("Mean profile"),
                               style = "padding: 20px;",
                               conditionalPanel(
                                 condition = "output.user_data_loaded_mean == false",
                                 div(
                                   class = "main-content",
                                   style = "text-align: center; padding: 100px 20px;",
                                   card(
                                     class = "unified-card",
                                     style = "max-width: 600px; margin: 0 auto;",
                                     card_header("No Mean Data Loaded"),
                                     card_body(
                                       icon("database", style = "font-size: 60px; color: #6c757d; margin-bottom: 20px;"),
                                       h4("Please upload your mean data first", style = "color: #6c757d;"),
                                       p("Go to the 'Plot' tab to load your BE data."),
                                       actionButton("go_to_upload",
                                                    label = div(icon("cloud-upload-alt"), "Go to Upload Mean Data"),
                                                    class = "unified-btn",
                                                    style = "margin-top: 20px;"
                                       )
                                     )
                                   )
                                 )
                               ),
                               conditionalPanel(
                                 condition = "output.user_data_loaded_mean == true",
                                 card(
                                   class = "unified-card card-accent-blue",
                                   card_header(span("Mean profile")),
                                   div(
                                     style = "display: flex; align-items: flex-end; gap: 15px; width: 100%;",
                                     
                                     # Внедряем кусочек CSS, который уберет скрытый отступ Bootstrap именно внутри этой группы
                                     tags$style(HTML("
    .no-margin-group .form-group { margin-bottom: 0 !important; }
    .no-margin-group .shiny-input-container { margin-bottom: 0 !important; }
  ")),
                                     
                                     div(
                                       class = "no-margin-group",
                                       style = "width: 30%;",
                                       radioGroupButtons(
                                         "y_coord_mean",
                                         "Select option:",
                                         choices = c("Normal" = "norm", "Log" = "ln"),
                                         status = "primary",
                                         justified = TRUE,
                                         width = "100%"
                                       )
                                     ),
                                     
                                     div(
                                       style = "width: 30%;",
                                       downloadButton(
                                         outputId = "download_mean_plot", 
                                         label = "Скачать средний профиль", 
                                         class = "btn-primary" ,
                                         style = "width: 100%; margin-bottom: 0 !important;"
                                       )
                                     )
                                   ),

                                   # Кнопка для скачивания усредненного профиля
                                   div(class = "unified-table", style = "padding: 10px;", plotOutput("mean_plot"))
                                 )
                               )
                             )
                           )
                         )
                       ))
    ), 
    tabPanel(
      title = 'Individual Concentrations',
      icon = icon("table"),
      conditionalPanel(
        condition = "output.user_data_loaded_mean == false",
        div(
          class = "main-content",
          style = "text-align: center; padding: 100px 20px;",
          card(
            class = "unified-card",
            style = "max-width: 600px; margin: 0 auto;",
            card_header("No Mean Data Loaded"),
            card_body(
              icon("database", style = "font-size: 60px; color: #6c757d; margin-bottom: 20px;"),
              h4("Please upload your mean data first", style = "color: #6c757d;"),
              p("Go to the 'Plot' tab to load your BE data."),
              actionButton("go_to_upload",
                           label = div(icon("cloud-upload-alt"), "Go to Upload Mean Data"),
                           class = "unified-btn",
                           style = "margin-top: 20px;"
              )
            )
          )
        )
      ),
      conditionalPanel(
        condition = "output.user_data_loaded_mean == true",
        div(style = "padding: 20px;", uiOutput("dynamic_concentrations_tab_ui")),
        
        card(
          class = "unified-card card-accent-blue",
          card_header(
            div(
              style = "display: flex; justify-content: space-between; align-items: center;",
              span("Individual concentrations")
            )
          ),
          
          div(class = "unified-table",
              style = "padding: 10px;",
              DTOutput('result_individ_conc')),
          
          div(class = "unified-table",
              style = "padding: 10px;",
              DTOutput('desc_statistic_individ_conc'))
        )
      )
      
    )
  )
)
## ==============================================================================## 3. СЕРВЕРНАЯ ЧАСТЬ (SERVER)## ==============================================================================
server <- function(input, output, session) {

  user_processed_data <- reactiveVal(NULL)
  user_processed_data_mean <- reactiveVal(NULL)
  user_processed_data_for_mean_profile <- reactiveVal(NULL)
  output$user_data_loaded <- reactive({ !is.null(user_processed_data()) })
  outputOptions(output, "user_data_loaded", suspendWhenHidden = FALSE)
  output$user_data_loaded_mean <- reactive({ !is.null(user_processed_data_mean()) })
  outputOptions(output, "user_data_loaded_mean", suspendWhenHidden = FALSE)
  output$user_data_for_mean_profile_loaded <- reactive({ !is.null(user_processed_data_for_mean_profile()) })
  outputOptions(output, "user_data_for_mean_profile_loaded", suspendWhenHidden = FALSE)
  data <- reactive({
    req(input$user_data)
    file_ext <- tools::file_ext(input$user_data$name)
    
    if (tolower(file_ext) == "csv") {
      df <- read.csv(input$user_data$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (tolower(file_ext) %in% c("xlsx", "xls")) {
      df <- read_excel(input$user_data$datapath)
    } else {
      return(NULL)
    }
    return(df)
  })

  output$data_rows <- renderText({ req(data()); paste(nrow(data()), "rows total") })
  observeEvent(input$load_user_data, {
    req(data())
    if (input$user_data$size > 50*1024^2) {
      output$upload_status <- renderUI({ div(class = "unified-alert unified-alert-danger", "File size error: Exceeds 50 MB.") })
      return()
    }
    user_processed_data(data())
    output$upload_status <- renderUI({ div(class = "unified-alert unified-alert-success", icon("check-circle"), tags$strong(" Data loaded successfully!")) })
  })
  ## --- Динамический инпут файла Mean на Вкладке 1 (Только для Replicate) ---


  data_mean <- reactive({
    req(input$user_data_mean)
    file_ext <- tools::file_ext(input$user_data_mean$name)
    
    if (tolower(file_ext) == "csv") {
      df <- read.csv(input$user_data_mean$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (tolower(file_ext) %in% c("xlsx", "xls")) {
      df <- read_excel(input$user_data_mean$datapath)
    } else {
      return(NULL)
    }
    return(df)
  })

  observeEvent(input$load_user_data_mean, {
    req(data_mean())
    user_processed_data_mean(data_mean())
    output$upload_status_mean <- renderUI({ div(class = "unified-alert unified-alert-success", icon("check-circle"), tags$strong(" Data for mean profile loaded successfully!")) })
  })
  output$data_preview <- renderDT({
    req(user_processed_data())
    datatable(head(user_processed_data() %>% mutate(Concentration = Concentration %>% as.numeric()), 10), options = list(pageLength = 10, dom = 't', ordering = FALSE), rownames = FALSE) %>% formatStyle(columns = seq_len(ncol(user_processed_data())), fontSize = '12px')
  })
  nca_data <- reactive({ req(user_processed_data()); nca_calculation(user_processed_data()) })
  output$resultNCA <- renderDT({
    req(nca_data(), input$tr)
    datatable(nca_data() %>% filter(Formulation == input$tr) %>% nca_rounding(input$timeround, input$pkround),
              extensions = 'Buttons', options = list(dom = 'Brt', buttons = c('copy', 'csv'), paging = FALSE, scrollY = "600px", scrollX = TRUE)) %>%
      formatRound(columns = c('TMAX','TLST','LAMZHL'), digits = input$timeround) %>%
      formatRound(columns = c('CMAX','CLST','CLSTP','AUCLST','AUCIFP'), digits = input$pkround) %>%
      formatRound(columns = 'LAMZ', digits = 4) %>% formatRound(columns = 'LAMZNPT', digits = 0)
  })
  output$nca_desc_statistic <- renderDT({
    req(nca_data(), input$tr)
    datatable(nca_data() %>% filter(Formulation == input$tr) %>% desc_statistic() %>% nca_rounding(input$timeround, input$pkround),
              extensions = 'Buttons', options = list(dom = 'Brt', buttons = c('copy', 'csv'), paging = FALSE, scrollY = "600px", scrollX = TRUE)) %>%
      formatRound(columns = c('TMAX','TLST','LAMZHL'), digits = input$timeround) %>%
      formatRound(columns = c('CMAX','CLST','CLSTP','AUCLST','AUCIFP'), digits = input$pkround) %>%
      formatRound(columns = 'LAMZ', digits = 4) %>% formatRound(columns = 'LAMZNPT', digits = 0)
  })
  ## --- ДИНАМИЧЕСКИЙ UI СЕЛЕКТОРОВ BE (РАЗЛИЧНЫЕ НАБОРЫ ДЛЯ РАЗНЫХ ФАЙЛОВ) ---
  output$dynamic_be_inputs <- renderUI({
    req(input$design_type)
    if (input$design_type == "replicate") {
      fluidRow(
        column(4, selectInput('pk', strong('Select PK'), choices = c('CMAX', 'AUCLST', 'AUCIFP'))),
        column(4, numericInput('alpha', strong('Input alpha value'), value = 0.05, min = 0, max = 1)),
        column(4, selectInput('design', strong('Select design'), choices = c('2x2x4', '2x2x3', '2x3x3')))
      )
    } else {
      fluidRow(
        column(12, selectInput('pk', strong('Select PK'), choices = c('CMAX', 'AUCLST', 'AUCIFP')))
      )
    }
  })
  be_data <- reactive({
    req(nca_data(), input$pk)
    if (input$design_type == "replicate") {
      req(input$alpha, input$design)
      be_result_replicate(input$pk, nca_data(), input$alpha, design = input$design)
    } else {
      be_result_cross(input$pk, nca_data())
    }
  })
  output$resultCI <- renderDT({
    req(be_data())
    datatable(be_data()$CIres, rownames = FALSE, extensions = 'Buttons', options = list(dom = 'Brt', buttons = c('copy', 'csv')))
  })
  output$anova_total <- renderDT({
    req(be_data())
    target_anova <- if(input$design_type == "replicate") be_data()$anova else be_data()$anova_total
    datatable(target_anova, extensions = 'Buttons', options = list(dom = 'Brt', buttons = c('copy', 'csv'))) %>%
      formatRound(columns = c('Sum Sq','Mean Sq','F value','Pr(>F)'), digits = 5) %>% formatRound(columns = c('Df'), digits = 0)
  })
  output$cv_power <- renderUI({
    req(be_data())
    if (input$design_type == "replicate") {
      tags$div(style = "padding: 20px; border-radius: 10px; background: #f8f9fa;",
               tags$span(style = "font-size: 20px; font-weight: bold; margin-right: 20px;", paste0('CVintra RR: ', round(be_data()$cv_r, 2), " %")),
               tags$span(style = "font-size: 20px; font-weight: bold; margin-right: 20px;", paste0('CVintra TT: ', round(be_data()$cv_t, 2), " %")),
               tags$span(style = "font-size: 20px; font-weight: bold;", paste0('Power: ', round(be_data()$power, 5))))
    } else {
      tags$div(style = "padding: 20px; border-radius: 10px; background: #f8f9fa;",
               tags$span(style = "font-size: 20px; font-weight: bold; margin-right: 20px;", paste0('CVintra: ', round(be_data()$cv, 2), " %")),
               tags$span(style = "font-size: 20px; font-weight: bold;", paste0('Power: ', round(be_data()$power, 5))))
    }
  })
  ## --- ДИНАМИЧЕСКИЙ ВВОД ПАРАМЕТРОВ САЙДБАРА ГРАФИКОВ ---

  calculated_plots <- reactive({
    req(user_processed_data(), input$y_coord)
    subjects <- unique(user_processed_data()$Subject)
    if (input$design_type == "replicate") {
      req(input$max_time, input$breaks)
      lapply(subjects, function(s) {
        subject_data <- user_processed_data() %>% filter(Subject == s)
        all_subject_plots_rep(subject_data, y_coord = input$y_coord, max_time = input$max_time, breaks = input$breaks)
      })
    } else {
      lapply(subjects, function(s) {
        all_subject_plots_cross(user_processed_data() %>% filter(Subject == s), y_coord = input$y_coord,  max_time = input$max_time, breaks = input$breaks)
      })
    }
  })
  observe({
    current_plots <- calculated_plots()
    for (i in seq_along(current_plots)) {
      local({
        my_i <- i
        plot_name <- paste0("plot_", my_i)
        output[[plot_name]] <- renderPlot({ calculated_plots()[[my_i]] })
      })
    }
  })
  output$plots_container <- renderUI({
    current_plots <- calculated_plots()
    plot_output_list <- lapply(seq_along(current_plots), function(i) {
      column(width = 6, plotOutput(paste0("plot_", i), height = "400px"), style = "margin-bottom: 20px;")
    })
    fluidRow(plot_output_list)
  })


  output$mean_plot <- renderPlot({
    req(user_processed_data_mean(), input$y_coord_mean)
    if (input$design_type == "replicate") {
      req(input$max_time, input$breaks)
      mean_subject_plots_rep(user_processed_data_mean(), y_coord = input$y_coord_mean, max_time = input$max_time, breaks = input$breaks)
    } else {
      req(input$max_time, input$breaks)
      mean_subject_plots_cross(user_processed_data_mean(), y_coord = input$y_coord_mean, max_time = input$max_time, breaks = input$breaks)
    }
  })
  output$overlap_plot <- renderPlot({
    req(user_processed_data(), input$y_coord_overlap, input$formulation)
    if (input$design_type == "replicate") {
      req(input$max_time, input$breaks)
      overlap_plot_rep(user_processed_data(), y_coord = input$y_coord_overlap, form = input$formulation, max_time = input$max_time, breaks = input$breaks)
    } else {
      req(input$max_time, input$breaks)
      overlap_plot_cross(user_processed_data(), y_coord = input$y_coord_overlap, form = input$formulation, max_time = input$max_time, breaks = input$breaks)
    }
  })
  
  # --- ДИНАМИЧЕСКИЙ КОНТЕНТ ВКЛАДКИ INDIVIDUAL CONCENTRATIONS (ТОЛЬКО ДЛЯ REPLICATE ИЛИ АДАПТИВНЫЙ ПОД ДАННЫЕ) ---
  output$dynamic_concentrations_tab_ui <- renderUI({
    req(input$design_type, user_processed_data_mean())
    target_df <- user_processed_data_mean()
    
    fluidRow(
      column(4, selectInput('prd', strong('Select period'), choices = target_df %>% pull(Period) %>% unique())),
      column(4, numericInput('conc_round', strong('Select concentration rounding'), value = 3, min = 0, max = 10))
    )


  })
  

  output$result_individ_conc <- renderDT({

    req(user_processed_data_mean(), input$prd)
    target_df <- user_processed_data_mean()

    data_individ <- target_df %>% filter(Period == input$prd) %>% data_individ_conc_calc() %>% individ_rounding(input$conc_round) %>% select(-Period)
    datatable(data_individ, extensions = 'Buttons', options = list(dom = 'Brt', buttons = c('copy', 'csv'), paging = FALSE, scrollY = "600px", scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(data_individ,is.numeric)), digits = input$conc_round)
  })
  
  output$desc_statistic_individ_conc <- renderDT({
    req(user_processed_data_mean(), input$prd)
    target_df <- user_processed_data_mean()

    data_individ_stat <- target_df %>% filter(Period == input$prd) %>% data_individ_conc_calc() %>% desc_statistic() %>% individ_rounding(input$conc_round)
    datatable(data_individ_stat, extensions = 'Buttons', options = list(dom = 'Brt', buttons = c('copy', 'csv'), paging = FALSE, scrollY = "600px", scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(data_individ_stat,is.numeric)), digits = input$conc_round)
  })
  ## --- Генерация примера файла ---

  
  # Скачивание примеров файлов crossover
  output$download_example_data_cross <- downloadHandler(
    filename = function() {
      "data_example.csv" # Возвращаем .csv, так как пишем через write.csv
    },
    content = function(file) {
      example_path <- "data/example/data_example.csv"
      if (file.exists(example_path)) {
        file.copy(example_path, file)
      } else {
        # Создаем реалистичный пример
        set.seed(123)
        n_subjects <- 20
        
        # Задаем сетку времени (до введения, частые точки для Cmax, редкие для элиминации)
        time_points <- c(-1.00, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 6.00, 8.00, 12.00, 14.00, 24.00)
        
        # Общее количество строк: субъект * период * временные точки
        n_periods <- 2 
        
        # Создаем каркас таблицы
        example_data <- expand.grid(
          Time = time_points,
          Period = 1:n_periods,
          Subject = sprintf("S%02d", 1:n_subjects)
        )
        
        # Назначаем Sequence и Formulation в зависимости от субъекта и периода (дизайн 2х2)
        example_data <- example_data %>%
          mutate(
            Subject_Num = as.numeric(gsub("S", "", Subject)),
            Sequence = ifelse(Subject_Num %% 2 == 0, "TR", "RT"),
            Formulation = case_when(
              Sequence == "TR" & Period == 1 ~ "T",
              Sequence == "TR" & Period == 2 ~ "R",
              Sequence == "RT" & Period == 1 ~ "R",
              Sequence == "RT" & Period == 2 ~ "T"
            )
          )
        
        # ГЕНЕРАЦИЯ ФК ПРОФИЛЕЙ (Одночастевая модель)
        example_data <- example_data %>%
          rowwise() %>%
          mutate(
            Ka = runif(1, 1.2, 1.8),    # Скорость всасывания
            Ke = runif(1, 0.1, 0.25),   # Скорость элиминации
            V  = runif(1, 40, 60),      # Объем распределения
            Dose = 500,                 # Условная доза
            
            F_bio = ifelse(Formulation == "T", runif(1, 0.85, 1.05), runif(1, 0.90, 1.10)),
            
            Concentration = if (Time <= 0) {
              0.000 
            } else {
              conc_val <- (Dose * F_bio * Ka) / (V * (Ka - Ke)) * (exp(-Ke * Time) - exp(-Ka * Time))
              pmax(0, conc_val + rnorm(1, mean = 0, sd = conc_val * 0.05)) %>% round(3)
            }
          ) %>%
          ungroup() %>%
          select(Subject, Time, Concentration, Formulation, Sequence, Period) %>%
          arrange(Subject, Period, Time)
        
        # Превращаем BLOQ некоторые нулевые значения
        example_data <- example_data %>%
          mutate(Concentration = ifelse(Time == -1.00, 0, as.character(Concentration)))
        
        write.csv(example_data, file, row.names = FALSE)
      }
    }
  )
  

  # Скачивание примеров файлов replicate
  output$download_example_data_rep <- downloadHandler(
    filename = function() {
      "data_example_replicate.csv" # Возвращаем .csv, так как пишем через write.csv
    },
    content = function(file) {
      example_path <- "data/example/data_example_replicate.csv"
      if (file.exists(example_path)) {
        file.copy(example_path, file)
      } else {
        # Создаем реалистичный пример
        set.seed(123)
        n_subjects <- 20
        
        # Задаем сетку времени (до введения, частые точки для Cmax, редкие для элиминации)
        time_points <- c(-1.00, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 6.00, 8.00, 12.00, 14.00, 24.00)
        
        # Общее количество строк: субъект * период * временные точки
        n_periods <- 4 
        
        # Создаем каркас таблицы
        example_data <- expand.grid(
          Time = time_points,
          Period = 1:n_periods,
          Subject = sprintf("S%02d", 1:n_subjects)
        )
        
        # Назначаем Sequence и Formulation в зависимости от субъекта и периода (дизайн 2х2)
        example_data <- example_data %>%
          mutate(
            Subject_Num = as.numeric(gsub("S", "", Subject)),
            Sequence = ifelse(Subject_Num %% 2 == 0, "TRTR", "RTRT"),
            Formulation = case_when(
              Sequence == "TRTR" & Period == 1 ~ "T",
              Sequence == "TRTR" & Period == 2 ~ "R",
              Sequence == "TRTR" & Period == 3 ~ "T",
              Sequence == "TRTR" & Period == 4 ~ "R",
              Sequence == "RTRT" & Period == 1 ~ "R",
              Sequence == "RTRT" & Period == 2 ~ "T",
              Sequence == "RTRT" & Period == 3 ~ "R",
              Sequence == "RTRT" & Period == 4 ~ "T"
            )
          )
        
        # ГЕНЕРАЦИЯ ФК ПРОФИЛЕЙ (Одночастевая модель)
        example_data <- example_data %>%
          rowwise() %>%
          mutate(
            Ka = runif(1, 1.2, 1.8),    # Скорость всасывания
            Ke = runif(1, 0.1, 0.25),   # Скорость элиминации
            V  = runif(1, 40, 60),      # Объем распределения
            Dose = 500,                 # Условная доза
            
            F_bio = ifelse(Formulation == "T", runif(1, 0.85, 1.05), runif(1, 0.90, 1.10)),
            
            Concentration = if (Time <= 0) {
              0.000 
            } else {
              conc_val <- (Dose * F_bio * Ka) / (V * (Ka - Ke)) * (exp(-Ke * Time) - exp(-Ka * Time))
              pmax(0, conc_val + rnorm(1, mean = 0, sd = conc_val * 0.05)) %>% round(3)
            }
          ) %>%
          ungroup() %>%
          select(Subject, Time, Concentration, Formulation, Sequence, Period) %>%
          arrange(Subject, Period, Time)
        
        # Превращаем BLOQ некоторые нулевые значения
        example_data <- example_data %>%
          mutate(Concentration = ifelse(Time == -1.00, 0, as.character(Concentration)))
        
        write.csv(example_data, file, row.names = FALSE)
      }
    }
  )
  
  # Скачивание индивидуальных графиков архивом
  output$download_all_plots <- downloadHandler(
    filename = function() {
      paste0("individual_plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # 1. Валидация данных
      plots_list <- req(calculated_plots())
      
      # 2. Создаем изолированную временную подпапку, чтобы файлы не перепутались
      tmp_sub_dir <- file.path(tempdir(), paste0("plots_", sample(1000:9999, 1)))
      dir.create(tmp_sub_dir, showWarnings = FALSE)
      
      files_to_zip <- c()
      
      # 3. Сохраняем каждый график
      for (i in seq_along(plots_list)) {
        plot_obj <- plots_list[[i]]
        
        # Безопасно вытаскиваем ID субъекта прямо из заголовка графика, 
        subj_id <- plot_obj$labels$title
        if (is.null(subj_id) || subj_id == "") subj_id <- paste0("index_", i)
        
        file_name <- paste0("Subject_", subj_id, ".png")
        full_path <- file.path(tmp_sub_dir, file_name)
        
        # Записываем файл на диск
        ggplot2::ggsave(full_path, plot = plot_obj, width = 7, height = 5, dpi = 300)
        files_to_zip <- c(files_to_zip, full_path)
      }
      
      # 4. Создаем ZIP-архив с помощью пакета zip (функция zipr)
      zip::zipr(zipfile = file, files = files_to_zip, recurse = FALSE, include_directories = FALSE)
      
      # Очищаем за собой временные файлы
      unlink(tmp_sub_dir, recursive = TRUE)
    }
  )


  # Скачивание усредненного профиля
  output$download_mean_plot <- downloadHandler(
    filename = function() {
      paste0("mean_profile_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(user_processed_data_for_mean_profile(), input$y_coord_mean, input$max_time, input$breaks)

      if (input$design_type == "replicate") {
        p <- mean_subject_plots_rep(user_processed_data_for_mean_profile(), y_coord = input$y_coord_mean, max_time = input$max_time, breaks = input$breaks)
      } else {
        p <- mean_subject_plots_cross(user_processed_data_for_mean_profile(), y_coord = input$y_coord_mean, max_time = input$max_time, breaks = input$breaks)
      }

      # Сохраняем в файл, который затребовал Shiny
      ggplot2::ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )


  # Скачивание профиля в наложении
  output$download_overlap_plot <- downloadHandler(
    filename = function() {
      paste0("overlap_profile_", input$formulation, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(data(), input$y_coord_overlap, input$formulation, input$max_time, input$breaks)

      if (input$design_type == "replicate") {
        p <- overlap_plot_rep(user_processed_data(), y_coord = input$y_coord_overlap, form = input$formulation, max_time = input$max_time, breaks = input$breaks)
      } else {
        p <- overlap_plot_cross(user_processed_data(), y_coord = input$y_coord_overlap, form = input$formulation, max_time = input$max_time, breaks = input$breaks)
      }


      ggplot2::ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )

  
  observeEvent(input$go_to_upload, { updateNavbarPage(session, "main_navbar", selected = "Upload Your Data") })
}
shinyApp(ui = ui, server = server)



  
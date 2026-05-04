setwd("/Users/saitoumarie/Library/CloudStorage/Dropbox/Norway/class/Bio322_2025/Bio322_analysis")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(scales)
})

##############
# Common setup
##############

cake <- c(
  "#3D1B1C",
  "#743118",
  "#BB6A40",
  "#F0D097",
  "#DC788C",
  "#E62E38"
)

border_col <- cake[1]
fill_col   <- cake[3]
title_size <- 12

theme_bio322_classic <- function(title_size = 12) {
  theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = title_size),
      legend.position = "bottom"
    )
}

theme_bio322_bw <- function(title_size = 12) {
  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = title_size),
      legend.position = "bottom"
    )
}

##############
# Data: resource access
##############
total_students_2025 <- 87
resources_agg <- read_csv("BIO322_resources_access.csv", show_col_types = FALSE)

resources_agg2 <- resources_agg %>%
  mutate(
    Resource = str_squish(Resource),
    Category = str_squish(Category),
    Total_students = total_students_2025,
    student_access_rate = Students / Total_students,
    views_per_accessing_student = Page_Views / Students,
    
    Category2 = case_when(
      Resource %in% c(
        "Course Home",
        "Course Assignments",
        "Course Announcements",
        "Course Quizzes",
        "Course Grades",
        "Course Syllabus"
      ) ~ "Course navigation",
      
      str_detect(
        Resource,
        regex(
          "Announcement|Schedule|Reminder|Deadline|Clarification|Q&A|Feedback|Appreciation|conflicts|finish early|course evaluation|Master|Upcoming|Small schedule change|Questions",
          ignore_case = TRUE
        )
      ) ~ "Course announcements",
      
      str_detect(Resource, regex("BIO322_2025_\\d|Module_3_lecture", ignore_case = TRUE)) ~
        "Lecture slides",
      
      str_detect(
        Resource,
        regex("Self-Learning Guide|Pre-course Information|BIO322_2025_guidance", ignore_case = TRUE)
      ) ~ "Preparatory guidance",
      
      str_detect(Resource, regex("Quiz", ignore_case = TRUE)) &
        str_detect(Resource, regex("Module 0|Module 1|Module 2|Module 3|Module 4|Final", ignore_case = TRUE)) ~
        "Quizzes",
      
      str_detect(
        Resource,
        regex(
          "Module4|Module 4|template|Final Presentation|Peer Review|grant application|Research plan|how to write",
          ignore_case = TRUE
        )
      ) ~ "Final project support",
      
      str_detect(
        Resource,
        regex(
          "Exercise|Assignment|submission|chromatogram|gRNA|Functional validation",
          ignore_case = TRUE
        )
      ) ~ "Exercises and submissions",
      
      Category == "Further Reading" ~ "Further reading",
      Category == "Model answers" ~ "Model answers",
      
      TRUE ~ "Other"
    ),
    
    Category2 = factor(
      Category2,
      levels = c(
        "Course navigation",
        "Course announcements",
        "Preparatory guidance",
        "Lecture slides",
        "Exercises and submissions",
        "Quizzes",
        "Final project support",
        "Further reading",
        "Model answers",
        "Other"
      )
    )
  )
##############
# Figure 1. Course evaluation scores
##############

evaluation <- read_csv("evaluation.csv", show_col_types = FALSE)

eval_long <- evaluation %>%
  pivot_longer(
    cols = -Item,
    names_to = "year_group",
    values_to = "value"
  ) %>%
  mutate(
    year  = str_extract(year_group, "^\\d{4}") %>% as.integer(),
    group = if_else(str_detect(year_group, "Faculty"), "Faculty", "BIO322")
  ) %>%
  select(Item, year, group, value) %>%
  pivot_wider(names_from = group, values_from = value) %>%
  mutate(delta = BIO322 - Faculty)

item_levels <- c(
  "Understanding of learning objectives",
  "Course structure and organization",
  "Contribution of lectures to learning",
  "Contribution of other learning activities",
  "Perceived learning outcome",
  "Overall satisfaction"
)

eval_long2 <- eval_long %>%
  mutate(Item = factor(Item, levels = item_levels))

band_2124 <- eval_long2 %>%
  filter(year >= 2021, year <= 2024) %>%
  group_by(Item) %>%
  summarise(
    ymin = min(BIO322, na.rm = TRUE),
    ymax = max(BIO322, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(ymin), is.finite(ymax))

pts_long <- eval_long2 %>%
  select(year, Item, BIO322, Faculty) %>%
  pivot_longer(
    cols = c(BIO322, Faculty),
    names_to = "Series",
    values_to = "Score"
  ) %>%
  mutate(
    Series = recode(Series, BIO322 = "BIO322", Faculty = "Faculty mean")
  )

shape_eval <- c("BIO322" = 21, "Faculty mean" = 24)
fill_eval  <- c("BIO322" = cake[2], "Faculty mean" = "#c1f4ca")

fig1_evaluation <- ggplot() +
  geom_rect(
    data = band_2124,
    aes(xmin = 2020.6, xmax = 2025.4, ymin = ymin, ymax = ymax),
    fill = "grey85"
  ) +
  geom_line(
    data = eval_long2,
    aes(x = year, y = BIO322),
    linewidth = 0.6,
    colour = border_col,
    show.legend = FALSE
  ) +
  geom_point(
    data = pts_long,
    aes(x = year, y = Score, shape = Series, fill = Series),
    size = 2.4,
    stroke = 1.0,
    colour = border_col,
    na.rm = TRUE
  ) +
  facet_wrap(~ Item, ncol = 3) +
  scale_x_continuous(breaks = 2021:2025) +
  scale_y_continuous(
    breaks = c(4, 5),
    limits = c(3.2, 5.5)
  ) +
  scale_shape_manual(values = shape_eval) +
  scale_fill_manual(values = fill_eval) +
  labs(
    x = "Year",
    y = "Mean score (1–6)",
    shape = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(hjust = 1),
    legend.position = "bottom"
  )

fig1_evaluation

ggsave("Figure1_course_evaluation.pdf", fig1_evaluation, width = 10, height = 4.5)
ggsave("Figure1_course_evaluation.png", fig1_evaluation, width = 10, height = 4.5, dpi = 300)

##############
# Figure 2. Preparatory materials before and after guide consolidation
##############

self_study <- read_csv("self_study.csv", show_col_types = FALSE)

shape_year <- c("2023" = 22, "2024" = 24, "2025" = 21)
fill_year  <- c("2023" = cake[4], "2024" = cake[3], "2025" = cake[5])

self_study2 <- self_study %>%
  mutate(
    Year = factor(Year, levels = c(2023, 2024, 2025)),
    is_guide = Content == "Self-Learning Guide",
    Content_lab = if_else(is_guide, "Guide", Content),
    content_num = suppressWarnings(as.numeric(Content)),
    x_ord = if_else(is_guide, -Inf, content_num),
    student_access_rate = Students / Total_students,
    views_per_accessing_student = Page_Views / Students
  )

self_study_levels <- self_study2 %>%
  distinct(Content_lab, x_ord) %>%
  arrange(x_ord) %>%
  pull(Content_lab)

self_study2 <- self_study2 %>%
  mutate(x = factor(Content_lab, levels = self_study_levels))

p_fig2_a <- ggplot(
  self_study2,
  aes(x = x, y = student_access_rate, shape = Year, fill = Year)
) +
  geom_point(
    aes(size = is_guide),
    stroke = 1.1,
    colour = border_col
  ) +
  scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 4.8), guide = "none") +
  scale_shape_manual(values = shape_year) +
  scale_fill_manual(values = fill_year) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x = NULL,
    y = "Student access rate",
    title = "Student access rate"
  ) +
  theme_bio322_classic(title_size)

p_fig2_b <- ggplot(
  self_study2,
  aes(x = x, y = views_per_accessing_student, shape = Year, fill = Year)
) +
  geom_point(
    aes(size = is_guide),
    stroke = 1.1,
    colour = border_col
  ) +
  scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 4.8), guide = "none") +
  scale_shape_manual(values = shape_year) +
  scale_fill_manual(values = fill_year) +
  scale_y_continuous(limits = c(0, 6), breaks = 0:6) +
  labs(
    x = NULL,
    y = "Page views per accessing student",
    title = "Views per accessing student"
  ) +
  theme_bio322_classic(title_size)

fig2_self_study <- (p_fig2_a + p_fig2_b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

fig2_self_study

ggsave("Figure2_self_learning_guide_access.pdf", fig2_self_study, width = 9, height = 3.8)
ggsave("Figure2_self_learning_guide_access.png", fig2_self_study, width = 9, height = 3.8, dpi = 300)

##############
# Figure 3. Resource access by category
##############

resources_plot <- resources_agg2 %>%
  filter(Category2 != "Other") %>%
  mutate(
    Category2_lab = case_when(
      Category2 == "Quizzes" ~ "Quizzes (mandatory/optional)",
      TRUE ~ as.character(Category2)
    )
  )

category_levels <- resources_plot %>%
  group_by(Category2_lab) %>%
  summarise(
    med_access_rate = median(student_access_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(med_access_rate)) %>%
  pull(Category2_lab)

category_levels_rev <- rev(category_levels)

resources_plot <- resources_plot %>%
  mutate(Category2_ord = factor(Category2_lab, levels = category_levels_rev))

resources_plot_left <- resources_plot

resources_plot_right <- resources_plot %>%
  filter(Resource != "Course Home")

p_fig3_a <- ggplot(
  resources_plot_left,
  aes(x = Category2_ord, y = student_access_rate)
) +
  geom_boxplot(outlier.shape = NA, fill = fill_col, colour = border_col) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  coord_flip() +
  theme_bio322_bw(title_size) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    y = "Student access rate",
    title = "Student access rate"
  )

p_fig3_b <- ggplot(
  resources_plot_right,
  aes(x = Category2_ord, y = views_per_accessing_student)
) +
  geom_boxplot(outlier.shape = NA, fill = fill_col, colour = border_col) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  coord_flip() +
  theme_bio322_bw(title_size) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    y = "Page views per accessing student",
    title = "Views per accessing student"
  )

fig3_category_overview <- (p_fig3_a + p_fig3_b) +
  plot_annotation(tag_levels = "A")

fig3_category_overview

ggsave(
  "Figure3_resource_category_overview.pdf",
  fig3_category_overview,
  width = 10,
  height = 4
)

ggsave(
  "Figure3_resource_category_overview.png",
  fig3_category_overview,
  width = 10,
  height = 4,
  dpi = 300
)

##############
# Figure 4. Lecture slide access across course sequence
##############
session_levels <- c("0.1", "1.1", "1.2", "1.3", "2.1", "2.2", "2.3", "3.1", "3.2", "3.3")

slides <- resources_agg2 %>%
  filter(Category2 == "Lecture slides") %>%
  mutate(
    session = str_extract(Resource, "\\d+\\.\\d+"),
    session = case_when(
      Resource == "Module_3_lecture_3_functional_genomics_2025.pdf" ~ "3.2",
      TRUE ~ session
    ),
    session = factor(session, levels = session_levels)
  ) %>%
  arrange(session)

p_fig4_a <- ggplot(slides, aes(x = session, y = student_access_rate, group = 1)) +
  geom_line(colour = border_col, linewidth = 0.5) +
  geom_point(size = 3.5, shape = 21, fill = fill_col, colour = border_col, stroke = 1.1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  theme_bio322_classic(title_size) +
  labs(
    x = NULL,
    y = "Student access rate",
    title = "Student access rate for lecture slides"
  )

p_fig4_b <- ggplot(slides, aes(x = session, y = views_per_accessing_student, group = 1)) +
  geom_line(colour = border_col, linewidth = 0.5) +
  geom_point(size = 3.5, shape = 21, fill = fill_col, colour = border_col, stroke = 1.1) +
  scale_y_continuous(limits = c(0, 6), breaks = 0:6) +
  theme_bio322_classic(title_size) +
  labs(
    x = NULL,
    y = "Page views per accessing student",
    title = "Views per accessing student"
  )

fig4_slides <- (p_fig4_a + p_fig4_b) +
  plot_annotation(tag_levels = "A")

fig4_slides


ggsave("Figure4_lecture_slides_access.pdf", fig4_slides, width = 8, height = 3.4)
ggsave("Figure4_lecture_slides_access.png", fig4_slides, width = 8, height = 3.4, dpi = 300)

##############
# Figure 5. Further-reading item access
##############
further_order <- tribble(
  ~fr_order, ~session, ~resource_key, ~reading_lab,
  1,  "0.1", "Next-generation sequencing technologies", "0.1a",
  2,  "0.1", "Genomics_biodiversity_conservation", "0.1b",
  3,  "1.1", "The importance of genomic variation", "1.1a",
  4,  "1.1", "Constraints and plasticity", "1.1b",
  5,  "1.2", "Ten things about transposable elements", "1.2a",
  6,  "1.2", "Taming transposable elements", "1.2b",
  7,  "1.3", "consequences_of_hybridization", "1.3a",
  8,  "2.1", "current population genomics methods", "2.1a",
  9,  "2.2", "Benefits and limitations of GWAS", "2.2a",
  10, "2.3", "A_new_era_in_functional_genomics_screens", "2.3a",
  11, "3.1", "RNA_sequencing_the_teenage_years", "3.1a",
  12, "3.2", "Chromosomal_instability", "3.2a",
  13, "3.2", "New insights into antibody structure", "3.2b",
  14, "3.3", "The relationship between genome structure and function", "3.3a",
  15, "4.0", "write a grant application", "4.0a"
)

reading_levels <- further_order %>%
  arrange(fr_order) %>%
  pull(reading_lab)

further_df <- further_order %>%
  rowwise() %>%
  mutate(
    matched = list(
      resources_agg2 %>%
        filter(str_detect(Resource, regex(resource_key, ignore_case = TRUE)))
    )
  ) %>%
  unnest(matched) %>%
  ungroup() %>%
  select(fr_order, session, reading_lab, Resource, Students, Page_Views, Total_students) %>%
  arrange(fr_order) %>%
  mutate(
    student_access_rate = Students / Total_students,
    views_per_accessing_student = Page_Views / Students,
    reading_lab = factor(reading_lab, levels = reading_levels)
  )

p_fig5_a <- ggplot(further_df, aes(x = reading_lab, y = student_access_rate)) +
  geom_col(fill = fill_col, colour = border_col, linewidth = 0.4) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bio322_classic(title_size) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    y = "Student access rate",
    title = "Student access rate for further-reading items"
  )

p_fig5_b <- ggplot(further_df, aes(x = reading_lab, y = views_per_accessing_student)) +
  geom_col(fill = fill_col, colour = border_col, linewidth = 0.4) +
  scale_y_continuous(
    breaks = 0:4
  ) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_bio322_classic(title_size) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    y = "Page views per accessing student",
    title = "Views per accessing student"
  )

fig5_further_reading <- (p_fig5_a + p_fig5_b) +
  plot_annotation(tag_levels = "A")

fig5_further_reading

ggsave("Figure5_further_reading_access.pdf", fig5_further_reading, width = 8, height = 4.2)
ggsave("Figure5_further_reading_access.png", fig5_further_reading, width = 8, height = 4.2, dpi = 300)
##############
# Figure 6. Quiz access
##############

quizzes <- resources_agg2 %>%
  filter(Category2 == "Quizzes") %>%
  mutate(
    quiz_order = case_when(
      str_detect(Resource, regex("Module 0", ignore_case = TRUE)) ~ 0,
      str_detect(Resource, regex("Module 1", ignore_case = TRUE)) ~ 1,
      str_detect(Resource, regex("Module 2", ignore_case = TRUE)) ~ 2,
      str_detect(Resource, regex("Module 3", ignore_case = TRUE)) ~ 3,
      str_detect(Resource, regex("Module 4", ignore_case = TRUE)) ~ 4,
      str_detect(Resource, regex("Final", ignore_case = TRUE)) ~ 5,
      TRUE ~ NA_real_
    ),
    quiz_type = case_when(
      quiz_order %in% c(0, 4) ~ "Optional",
      quiz_order %in% c(1, 2, 3, 5) ~ "Mandatory",
      TRUE ~ NA_character_
    ),
    Resource_lab = case_when(
      quiz_order == 0 ~ "0",
      quiz_order == 1 ~ "1",
      quiz_order == 2 ~ "2",
      quiz_order == 3 ~ "3",
      quiz_order == 4 ~ "4",
      quiz_order == 5 ~ "Final",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(quiz_order)) %>%
  mutate(
    Resource_lab = factor(
      Resource_lab,
      levels = c("0", "1", "2", "3", "4", "Final")
    ),
    quiz_type = factor(quiz_type, levels = c("Mandatory", "Optional")),
    student_access_rate = Students / Total_students,
    views_per_accessing_student = Page_Views / Students
  ) %>%
  arrange(quiz_order)

quiz_fill <- c(
  "Mandatory" = cake[3],
  "Optional" = "#F0D097"
)

p_fig6_a <- ggplot(
  quizzes,
  aes(x = Resource_lab, y = student_access_rate, fill = quiz_type)
) +
  geom_col(colour = border_col, linewidth = 0.4) +
  scale_fill_manual(values = quiz_fill) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  theme_bio322_classic(title_size) +
  labs(
    x = NULL,
    y = "Student access rate",
    fill = NULL,
    title = "Student access rate for quizzes"
  )

p_fig6_b <- ggplot(
  quizzes,
  aes(x = Resource_lab, y = views_per_accessing_student, fill = quiz_type)
) +
  geom_col(colour = border_col, linewidth = 0.4) +
  scale_fill_manual(values = quiz_fill) +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = 0:10
  ) +
  theme_bio322_classic(title_size) +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Page views per accessing student",
    fill = NULL,
    title = "Views per accessing student"
  )

fig6_quizzes <- (p_fig6_a + p_fig6_b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

fig6_quizzes

ggsave(
  "Figure6_quiz_access.pdf",
  fig6_quizzes,
  width = 8,
  height = 3.8
)

ggsave(
  "Figure6_quiz_access.png",
  fig6_quizzes,
  width = 8,
  height = 3.8,
  dpi = 300
)





###
##############
# Descriptive statistics for Figures 1–6
##############

total_students_2025 <- 87

##############
# Figure 1: Course evaluation
##############

fig1_stats <- eval_long2 %>%
  group_by(Item) %>%
  summarise(
    bio322_2025 = BIO322[year == 2025],
    pre2025_min = min(BIO322[year >= 2021 & year <= 2024], na.rm = TRUE),
    pre2025_median = median(BIO322[year >= 2021 & year <= 2024], na.rm = TRUE),
    pre2025_max = max(BIO322[year >= 2021 & year <= 2024], na.rm = TRUE),
    faculty_2025 = Faculty[year == 2025],
    above_pre2025_range = bio322_2025 > pre2025_max,
    delta_2025_vs_pre2025_median = bio322_2025 - pre2025_median,
    delta_2025_vs_faculty = bio322_2025 - faculty_2025,
    .groups = "drop"
  ) %>%
  mutate(
    figure = "Figure 1",
    resource_type = "Course evaluation"
  ) %>%
  select(
    figure,
    resource_type,
    Item,
    bio322_2025,
    pre2025_min,
    pre2025_median,
    pre2025_max,
    above_pre2025_range,
    delta_2025_vs_pre2025_median,
    faculty_2025,
    delta_2025_vs_faculty
  )

fig1_summary <- fig1_stats %>%
  summarise(
    figure = "Figure 1",
    resource_type = "Course evaluation",
    n_items = n(),
    n_above_pre2025_range = sum(above_pre2025_range, na.rm = TRUE),
    median_delta_2025_vs_pre2025_median = median(delta_2025_vs_pre2025_median, na.rm = TRUE),
    min_delta_2025_vs_pre2025_median = min(delta_2025_vs_pre2025_median, na.rm = TRUE),
    max_delta_2025_vs_pre2025_median = max(delta_2025_vs_pre2025_median, na.rm = TRUE),
    median_delta_2025_vs_faculty = median(delta_2025_vs_faculty, na.rm = TRUE),
    min_delta_2025_vs_faculty = min(delta_2025_vs_faculty, na.rm = TRUE),
    max_delta_2025_vs_faculty = max(delta_2025_vs_faculty, na.rm = TRUE)
  )

##############
# Helper function for resource-access figures
##############

summarise_access <- function(df, figure_name, resource_type, total_students = total_students_2025) {
  df %>%
    mutate(
      student_access_rate = Students / total_students,
      views_per_accessing_student = Page_Views / Students
    ) %>%
    summarise(
      figure = figure_name,
      resource_type = resource_type,
      n_items = n(),
      
      min_access_rate = min(student_access_rate, na.rm = TRUE),
      median_access_rate = median(student_access_rate, na.rm = TRUE),
      max_access_rate = max(student_access_rate, na.rm = TRUE),
      n_access_rate_above_50pct = sum(student_access_rate > 0.5, na.rm = TRUE),
      
      min_views_per_accessing_student = min(views_per_accessing_student, na.rm = TRUE),
      median_views_per_accessing_student = median(views_per_accessing_student, na.rm = TRUE),
      max_views_per_accessing_student = max(views_per_accessing_student, na.rm = TRUE)
    )
}

##############
# Figure 2: Preparatory materials
##############

fig2_stats <- self_study2 %>%
  mutate(
    student_access_rate = Students / Total_students,
    views_per_accessing_student = Page_Views / Students,
    resource_group = case_when(
      is_guide ~ "2025 self-learning guide",
      Year %in% c("2023", "2024") ~ "2023–2024 session-linked preparatory items",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(resource_group) %>%
  summarise(
    figure = "Figure 2",
    resource_type = first(resource_group),
    n_items = n(),
    
    min_access_rate = min(student_access_rate, na.rm = TRUE),
    median_access_rate = median(student_access_rate, na.rm = TRUE),
    max_access_rate = max(student_access_rate, na.rm = TRUE),
    n_access_rate_above_50pct = sum(student_access_rate > 0.5, na.rm = TRUE),
    
    min_views_per_accessing_student = min(views_per_accessing_student, na.rm = TRUE),
    median_views_per_accessing_student = median(views_per_accessing_student, na.rm = TRUE),
    max_views_per_accessing_student = max(views_per_accessing_student, na.rm = TRUE),
    .groups = "drop"
  )

##############
# Figure 3: Resource categories
##############

fig3_stats <- resources_plot %>%
  mutate(
    student_access_rate = Students / total_students_2025,
    views_per_accessing_student = Page_Views / Students
  ) %>%
  filter(Resource != "Course Home") %>%
  group_by(Category2_lab) %>%
  summarise(
    figure = "Figure 3",
    resource_type = first(Category2_lab),
    n_items = n(),
    
    min_access_rate = min(student_access_rate, na.rm = TRUE),
    median_access_rate = median(student_access_rate, na.rm = TRUE),
    max_access_rate = max(student_access_rate, na.rm = TRUE),
    n_access_rate_above_50pct = sum(student_access_rate > 0.5, na.rm = TRUE),
    
    min_views_per_accessing_student = min(views_per_accessing_student, na.rm = TRUE),
    median_views_per_accessing_student = median(views_per_accessing_student, na.rm = TRUE),
    max_views_per_accessing_student = max(views_per_accessing_student, na.rm = TRUE),
    .groups = "drop"
  )

course_home_stat <- resources_plot %>%
  filter(Resource == "Course Home") %>%
  transmute(
    figure = "Figure 3",
    resource_type = "Course Home",
    Students,
    Page_Views,
    student_access_rate = Students / total_students_2025,
    views_per_accessing_student = Page_Views / Students
  )

##############
# Figure 4: Lecture slides
##############

fig4_stats <- summarise_access(
  df = slides,
  figure_name = "Figure 4",
  resource_type = "Lecture slides"
)

##############
# Figure 5: Further reading
##############

fig5_stats_all <- summarise_access(
  df = further_df,
  figure_name = "Figure 5",
  resource_type = "Further reading, all items"
)

fig5_stats_no_4 <- summarise_access(
  df = further_df %>% filter(session != "5"),
  figure_name = "Figure 5",
  resource_type = "Further reading, excluding 5.0 final-project item"
)

fig5_stats_4 <- summarise_access(
  df = further_df %>% filter(session == "5"),
  figure_name = "Figure 5",
  resource_type = "Further reading, 5.0 final-project item"
)

##############
# Figure 6: Quizzes
##############

fig6_stats_all <- summarise_access(
  df = quizzes,
  figure_name = "Figure 6",
  resource_type = "Quizzes, all"
)

fig6_stats_by_type <- quizzes %>%
  mutate(
    student_access_rate = Students / total_students_2025,
    views_per_accessing_student = Page_Views / Students
  ) %>%
  group_by(quiz_type) %>%
  summarise(
    figure = "Figure 6",
    resource_type = paste0("Quizzes, ", first(quiz_type)),
    n_items = n(),
    
    min_access_rate = min(student_access_rate, na.rm = TRUE),
    median_access_rate = median(student_access_rate, na.rm = TRUE),
    max_access_rate = max(student_access_rate, na.rm = TRUE),
    n_access_rate_above_50pct = sum(student_access_rate > 0.5, na.rm = TRUE),
    
    min_views_per_accessing_student = min(views_per_accessing_student, na.rm = TRUE),
    median_views_per_accessing_student = median(views_per_accessing_student, na.rm = TRUE),
    max_views_per_accessing_student = max(views_per_accessing_student, na.rm = TRUE),
    .groups = "drop"
  )

##############
# Combined summary tables
##############

resource_access_summary <- bind_rows(
  fig2_stats,
  fig3_stats,
  fig4_stats,
  fig5_stats_all,
  fig5_stats_no_4,
  fig5_stats_4,
  fig6_stats_all,
  fig6_stats_by_type
) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )

resource_access_summary

fig1_summary

fig1_stats

course_home_stat

write_csv(fig1_stats, "Figure1_course_evaluation_item_stats.csv")
write_csv(fig1_summary, "Figure1_course_evaluation_summary.csv")
write_csv(resource_access_summary, "Figures2_to_6_resource_access_summary.csv")
write_csv(course_home_stat, "Figure3_course_home_stat.csv")
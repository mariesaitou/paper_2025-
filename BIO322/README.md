BIO322 Course Redesign Analysis
# Overview

This repository contains the analysis code and summary outputs for the BIO322 graduate genomics course redesign study.

The analysis focuses on course-level evaluation of a redesign that aimed to improve course navigability by centralizing and standardizing preparatory materials, without modifying core content or assessment structure.

The results are based on:

Institutional course evaluation data (2021–2025)
Aggregated learning management system (LMS) resource-access data (2023–2025)

No individual-level data are included.

# Data

The analysis requires the following input files (not included):

evaluation.csv
Course evaluation scores across years and items
BIO322_resources_access.csv
Aggregated LMS resource access data (students and page views)
self_study.csv
Access data for preparatory materials (before and after guide consolidation)

All datasets are anonymized and aggregated. No identifiable student information is included.

Analysis structure

The script performs the following steps:

1. Data preprocessing
Cleaning and standardization of resource names and categories
Calculation of:
Student access rate
Views per accessing student
Recategorization of LMS resources into analysis-relevant groups
2. Figures

The script generates the following figures:

Figure 1: Course evaluation scores (BIO322 vs faculty mean, 2021–2025)
Figure 2: Access to preparatory materials before and after guide consolidation
Figure 3: Resource access by category
Figure 4: Lecture slide access across course sequence
Figure 5: Further-reading item access
Figure 6: Quiz access (mandatory vs optional)

All figures are saved as both .pdf and .png.

3. Descriptive statistics

The script outputs summary statistics corresponding to Figures 1–6, including:

Access rate distributions
Views per accessing student
Changes across years (Figure 1)
Category-level summaries

Exported files:

Figure1_course_evaluation_item_stats.csv
Figure1_course_evaluation_summary.csv
Figures2_to_6_resource_access_summary.csv
Figure3_course_home_stat.csv
Requirements

R (≥ 4.0) with the following packages:

tidyverse
readr
lubridate
stringr
patchwork
scales



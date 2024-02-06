### set the environment
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(gtsummary)
library(flextable)

### load data
data(CO2)
head(CO2)

### summary table
tbl <- CO2 %>% 
  select(!Plant) %>% 
  tbl_summary(by=Type,
              type = list(conc ~ "continuous"),
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} / {N} ({p}%)"),
              digits = all_continuous() ~ 2,
              label = list(Treatment~ "Hello World")) %>% 
  add_p() %>% 
  add_overall() %>%
  modify_spanning_header(c("stat_1","stat_2") ~ "**Location**") %>%
  modify_header(label="**Variables**") %>%
  modify_footnote(label ~ "variables", stat_0 ~ NA) %>%
  as_flex_table() %>%
  set_caption(as_paragraph(as_b("Table 1.1 CO2")))

  # modify_table_body()
  # modify_column_indent()
  
save_as_docx(tbl, path = 'summary_table.docx')

# table 2
tbl2 <- CO2 %>% 
  tbl_cross(row = Type,
            col = Treatment,
            percent = "cell") %>% 
  add_p() %>%
  as_flex_table()

save_as_docx(tbl2, path = 'summary_table2.docx')

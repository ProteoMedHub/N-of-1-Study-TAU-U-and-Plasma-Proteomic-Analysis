# title: "N-of-1 HOMDD-1 Studies TAU-U Analysis"
author: Antonio Brazil Viana Junior 
format: html
editor: visual
---

```{r}
#| echo: true
#| code-fold: true
#| code-summary: "expand for full code"
# IMPORT ------------------------------------------------------------------
library(readxl)
data <- read_excel("data.xlsx") #Replace the file name “data.xlsx” with the name of your Excel file.

pacman::p_load(tidyverse,scan,kableExtra,SingleCaseES)
```

### Plots

```{r}
#| echo: true
#| code-fold: true
#| code-summary: "expand for full code"

data %>% filter(!is.na(`Treatment Period`)) %>% 
  group_by(`Treatment Period`) %>% 
  mutate(n=row_number()) %>% 
  ungroup() %>% 
  ggplot( aes(n, PCS)) + 
  geom_line(size=0.5) +
  geom_point()+
  geom_smooth(method = "lm",se = F,linetype = "dashed", color="grey")+
  xlab("\nSampling") + ylab("PCS\n") +
  theme_classic() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) + facet_grid(~`Treatment Period`)

data %>% filter(!is.na(`Treatment Period`)) %>% 
  group_by(`Treatment Period`) %>% 
  mutate(n=row_number()) %>% 
  ungroup() %>% 
  ggplot( aes(n, MCS)) + 
  geom_line(size=0.5) +
  geom_point()+
  geom_smooth(method = "lm",se = F,linetype = "dashed", color="gray")+
  xlab("\nSampling") + ylab("MCS\n") +
  theme_classic() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) + facet_grid(~`Treatment Period`)

data %>% filter(!is.na(`Treatment Period`)) %>% 
  group_by(`Treatment Period`) %>% 
  mutate(n=row_number()) %>% 
  ungroup() %>% 
  ggplot( aes(n, BDI)) + 
  geom_line(size=0.5) +
  geom_point()+
  geom_smooth(method = "lm",se = F,linetype = "dashed", color="gray")+
  xlab("\nSampling") + ylab("BDI\n") +
  theme_classic() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) + facet_grid(~`Treatment Period`)
```

## TAU-U

```{r}
#| echo: true
#| code-fold: true
#| code-summary: "expand for full code"
#Adjusting the parameters
data %>% arrange(`Treatment Period`)
n_a <- as.numeric(count(banco[data$treatment period =="A",])) 
n_b <- as.numeric(count(banco[data$treatment period =="B",]))
BD.Ordenado <- banco %>% arrange(treatment period)
vPCS <- as.vector(BD.Ordenado$PCS)
vBDI <- as.vector(BD.Ordenado$BDI)
vMCS <- as.vector(BD.Ordenado$MCS)


tPCS <- scdf(vPCS, phase_design = c(A = n_a, B = n_b), name = "PCS")
Result.tPCS <- tau_u(tPCS, method = "parker", tau_method = "b",ci = 0.95,ci_method="s") 
print(Result.tPCS, complete = TRUE) 


tBDI <- scdf(vBDI, phase_design = c(A = n_a, B = n_b), name = "BDI")
Result.tBDI <- tau_u(tBDI, method = "parker", tau_method = "b",ci = 0.95,ci_method="s") 
print(Result.tBDI, complete = TRUE)


tMCS <- scdf(vMCS, phase_design = c(A = n_a, B = n_b), name = "MCS")
Result.tMCS <- tau_u(tMCS, method = "parker", tau_method = "b",ci = 0.95,ci_method="s") 
print(Result.tMCS, complete = TRUE)
```

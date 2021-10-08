---
title: "Cnidarian Gene Expression Literature Review"
author: "Carly Kenkel & Rachel Wright"
date: "8/16/2021"
output:
  html_document: default
  pdf_document: default
---

# Load packages
```{r load packages, warning=F}
library(tidyverse) # for data wrangling, ggplots
library(cowplot) # for combining/modifying ggplots
library(networkD3) # for interactive Sankey plots
library(ggsci) # for color palettes
library(plotly) # for making interactive plots
library(UpSetR) # for making sunburst plots
library(orca) # for saving static plots from plotly
```

# Load literature review data sheet

```{r load lit review spreadsheet}
lit0 <- read.csv("LitReview.xlsx - primarLit.csv")
nrow(lit0)
```

Get rid of non quantitative gene expression methods

```{r only GE}
levels(as.factor(lit0$Method))

lit0 <- lit0 %>% filter(!Method %in% c("Expressed sequence tags", "protein immunoblotting", 
                                     "reannealing rate assay", "Northern blot", "in situ hybridization"))
levels(as.factor(lit0$Method))
```
How many publications now?

```{r}
nrow(lit0)
```

# Make a separate column for experimental factors that are stress-related
Tally the number of stressors

```{r}
levels(lit0$Expt.Factors)
nonstress_factors <- list("developmental stage", "color morph", "source population", "time course", "allorecognition",
                       "antibiotics", "cryopreservation", "depth", "field sites", "genotype", 
                       "reproductive cycle", "settlement cue", "symbiosis", "symbiotic state",
                       "light:dark", "strain", "acclimation", "species", "tissue", "lunar cycle", 
                       "algal strain", "life stage", "feeding", "algal species")
nonstress_factors <- paste(unlist(nonstress_factors), collapse = "|")

lit0$Expt.Factors.Stress <- str_remove_all(lit0$Expt.Factors, nonstress_factors)
lit0$Expt.Factors.Stress <- str_remove(lit0$Expt.Factors.Stress, "^;+")
lit0$Expt.Factors.Stress <- str_remove(lit0$Expt.Factors.Stress, ";+$")
lit0$Expt.Factors.Stress[lit0$Expt.Factors.Stress==""] <- NA
levels(as.factor(lit0$Expt.Factors.Stress))
lit0$NumStress <- str_count(lit0$Expt.Factors.Stress, ";")+1
```

Rename some factors to match (e.g., temp and acidification is the same as acidification and temp)

```{r}
lit0 %>% select(Expt.Factors.Stress) %>% filter(grepl("acidification",lit0$Expt.Factors.Stress)) %>% distinct

lit0$Expt.Factors.Stress <- gsub("acidification;temperature","temperature;acidification",lit0$Expt.Factors.Stress)

lit0 %>% select(Expt.Factors.Stress) %>% filter(grepl("acidification",lit0$Expt.Factors.Stress)) %>% distinct

```

The interactions are always between two stressors.

Sometimes they are fully factorial and get counted as "3" stresses (e.g., temperature;salinity;temp x salinity).

Sometimes they are not fully factorial and only considered the 2 stresses together (e.g., temp x salinity).

For now I'm going to count either case as 2 stressors because that is the number of experimental variables that were manipulated.

```{r}
lit0$Interaction <- ifelse(grepl(" x ", lit0$Expt.Factors.Stress), 
                           "interaction", "none")
lit0$NumStressInt <- if_else(lit0$Interaction=="interaction",
                            true = 2, false = lit0$NumStress)

```

# Tally the number of experimental factors to distinguish between single- and multiple-stressor studies

```{r}
lit0  %>% filter(NumStressInt == 1) %>% summarize(n = n(),
                                                 prop = n()/nrow(lit0))
```

Number of multiple stressor studies

```{r}
lit0 %>% filter(NumStressInt > 1) %>% summarize(n = n(),
                                                 prop = n()/nrow(lit0))
```

Histogram of number of stressors

```{r}
lit0 %>%
  ggplot(aes(x = NumStressInt)) +
  geom_histogram(binwidth = 1, color="white") +
  xlab("Number of Stressors Examined")+
  ylab("Number of Publications")+
  theme_classic()
```

What percentage of studies include an interaction between multiple simultaneous stressors?

```{r}
lit0 %>% filter(!is.na(Expt.Factors.Stress)) %>% tally()
lit0 %>% filter(!is.na(Expt.Factors.Stress)) %>% filter(Interaction=="interaction") %>% tally()
```

What factors do those studies that include an interaction consider?

```{r}
lit0 %>% filter(Interaction=="interaction") %>% select(Expt.Factors.Stress)
```

# What proportion of metagenomic?

```{r}
lit0 %>% select(Scope) %>% group_by(Scope) %>% summarize(n())
```

# Format literature spreadsheet to save as a supplemental table

```{r}
# write.csv(lit0, file = "TableS1.csv", row.names = F)
```


# Format literature spreadsheet for analysis

Keep columns containing title (TI), Journal (SO), Scope, Method, DevelopmentalStage, Study.Focus, Expt.Factors, Expt.Factors.Stress, NumStressInt, Conservation.Relevance, Partner, FocalSpp, Focal.Fam, and year (PY)

```{r select columns}
lit1 <- lit0 %>% select(Scope, Method, DevelopmentalStage, Study.Focus, Expt.Factors, Expt.Factors.Stress, NumStressInt, 
                       Conservation.Relevance, Partner, FocalSpp, Focal.Fam,
                       PY, TI, SO)
```

How many total citations?

```{r tally citations unfiltered}
nrow(lit1)
```

Make new rows for each entry within the lit review terms

```{r individual rows}
lit1 <- lit1 %>% separate_rows(Method, sep = ";") %>% separate_rows(Study.Focus, sep = ";") %>%
  separate_rows(Expt.Factors, sep = ";") %>% separate_rows(Conservation.Relevance, sep = ";") %>%
  separate_rows(DevelopmentalStage, sep = ";")
```

Make new rows for each entry for a particular species. 
Match the species to the family

```{r}
lit1 <- lit1 %>% separate_rows(FocalSpp, sep = ";")
# write.csv(file="species.csv",levels(as.factor(lit$FocalSpp)))
# I manually added family taxa to the species we have here
spp2fam <- read.csv("species.csv")
spp2fam
levels(as.factor(spp2fam$family))

# test it out...
test <- merge(lit1, spp2fam, by.x = "FocalSpp", by.y = "species", all.x=T) %>%
  select(TI, DevelopmentalStage, FocalSpp, Focal.Fam, family)
test %>% select(TI) %>% distinct() %>% tally()

lit <-  merge(lit1, spp2fam, by.x = "FocalSpp", by.y = "species", all.x=T)
```

How many studies now?

```{r}
lit %>% select(TI) %>% distinct() %>% tally()
```

Make sure that every entry is named exactly the same (e.g., long term stress and long term stress (>31d) need to be renamed to be the same)

```{r check levels}
levels(lit$Scope)
levels(lit$DevelopmentalStage)
levels(as.factor(lit$Method))
levels(as.factor(lit$Study.Focus))
levels(as.factor(lit$Expt.Factors))
levels(as.factor(lit$Conservation.Relevance))
levels(as.factor(lit$Partner))
levels(as.factor(lit$FocalSpp))
levels(as.factor(lit$family))
```

Combine all "Symbiodinium"-type partners

```{r rename Symbiodinium}
levels(as.factor(lit$Partner))
lit$Partner <- str_replace(lit$Partner, "Symbiodinium microadriaticum", "Symbiodiniaceae")
lit$Partner <- str_replace(lit$Partner, "Symbiodinium", "Symbiodiniaceae")
```

# Save/Load: You can start here with the "lit" object

```{r START HERE}
# save(lit, file="litSearch.RData")
load("litSearch.RData")
```

# Method & Partner

Make a bar chart of study method

```{r bar chart study method}
lit %>% group_by(Method) %>% filter(!is.na(Method)) %>% distinct(TI) %>% tally() %>%
 ggplot() +
  geom_bar(aes(reorder(Method,n),n), stat="identity", width=1, color="white") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("")+
  ylab("Number of Studies")
```

Plot method over time

```{r plot method by year line, warning=F}
colors <- pal_npg()(6)
lit %>% group_by(Method, PY) %>% filter(!is.na(Method)) %>% distinct(TI) %>% tally() %>%
  ggplot() +
  geom_line(aes(x=PY,y=n, group=Method,color=Method), size=1.5) +
  scale_color_manual(values=colors)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("")+
  ylab("Number of Studies")
```

```{r plot method by year bar, warning=F}
colors <- pal_npg()(6)
lit %>% group_by(Method, PY) %>% filter(!is.na(Method)) %>% distinct(TI) %>% tally() %>%
  ggplot() +
  geom_bar(aes(x=PY,y=n, group=Method,fill=Method),position="stack", stat="identity") +
  theme_classic()+
  scale_fill_manual("Method",values=colors)+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.text=element_text(size=rel(0.75)),
        legend.title=element_text(size=rel(0.75)))+
  xlab("")+
  ylab("Number of Studies")
```

Remove SSH

```{r plot method by year bar remove SSH, warning=F}
colors <- pal_npg()(5)
plot.method.time <- lit %>% group_by(Method, PY) %>% 
  filter(!is.na(Method)) %>% filter(!Method=="SSH") %>% 
  distinct(TI) %>% tally() %>%
  ggplot() +
  geom_bar(aes(x=PY,y=n, group=Method,fill=Method),position="stack", stat="identity") +
  theme_classic()+
  scale_fill_manual("Method",values=colors)+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.text=element_text(size=rel(0.75)),
        legend.title=element_text(size=rel(0.75)))+
  xlab("")+
  ylab("Number of Studies")
plot.method.time
```

```{r plot partner by year, warning=F}
colors <- pal_npg()(9)
lit %>% group_by(Partner, PY) %>% filter(!is.na(Partner)) %>% distinct(TI) %>% tally() %>%
  ggplot() +
  geom_line(aes(x=PY,y=n, group=Partner,color=Partner), size=1.5) +
  theme_classic()+
  scale_color_manual(values=colors)+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.text=element_text(size=rel(0.75)),
        legend.title=element_text(size=rel(0.75)))+
  xlab("")+
  ylab("Number of Studies")
```

```{r plot partner by year bar remove rare methods, warning=F}
my_levels <- c("Coral host", 
                "Anemone host",
                "Coral host;Symbiodiniaceae", 
                "Coral host;Symbiodiniaceae;other eukaryotes",
                "Anemone host;Symbiodiniaceae",
                "Symbiodiniaceae",
                "Bacteria")

colors <- c("#E64B35FF", "#F39B7FFF", "#3C5488FF", "#8491B4FF",
            "#4DBBD5FF", "#00A087FF", "#91D1C2FF")

plot.partner.time <- lit %>% group_by(Partner) %>% filter(n()>5) %>% ungroup() %>%
  group_by(Partner, PY) %>% filter(!is.na(Partner)) %>% distinct(TI) %>% tally() %>%
  ggplot() +
  geom_bar(aes(x=PY,y=n, 
               group=factor(Partner,levels=my_levels),
               fill=factor(Partner,levels=my_levels)), 
           position="stack", stat="identity") +
  theme_classic()+
  scale_fill_manual(values=colors)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.text=element_text(size=rel(0.75)),
        legend.title=element_text(size=rel(0.75)))+  
  xlab("")+
  ylab("Number of Studies")
plot.partner.time
```

# Conservation Relevance

```{r plot conservation relevance, warning=F}
colors <- pal_npg()(10)
plot.conservation.time <- lit %>% group_by(Conservation.Relevance) %>% filter(n()>2) %>% ungroup() %>%
  group_by(Conservation.Relevance, PY) %>% filter(!is.na(Conservation.Relevance)) %>% distinct(TI) %>% tally() %>%
  ggplot() +
  geom_bar(aes(x=PY,y=n, group=Conservation.Relevance,
               fill=Conservation.Relevance), position="stack", stat="identity") +
  theme_classic()+
  scale_fill_manual("Conservation Relevance", values=colors)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.text=element_text(size=rel(0.75)),
        legend.title=element_text(size=rel(0.75)))+
  xlab("")+
  ylab("Number of Studies")
plot.conservation.time
```

Remove "management decisions" (vague)

```{r plot conservation relevance 2, warning=F}
colors <- pal_npg()(10)
plot.conservation.time <- lit %>% filter(!Conservation.Relevance=="management decisions") %>%
  group_by(Conservation.Relevance) %>% filter(n()>2) %>% ungroup() %>%
  group_by(Conservation.Relevance, PY) %>% filter(!is.na(Conservation.Relevance)) %>% distinct(TI) %>% tally() %>%
  ggplot() +
  geom_bar(aes(x=PY,y=n, group=Conservation.Relevance,
               fill=Conservation.Relevance), position="stack", stat="identity") +
  theme_classic()+
  scale_fill_manual("Conservation Relevance", values=colors)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.text=element_text(size=rel(0.75)),
        legend.title=element_text(size=rel(0.75)))+
  xlab("")+
  ylab("Number of Studies")
plot.conservation.time
```


Plot method, partner, and conservation relevance together

```{r barplot method and partner and conservation}
fig_stackedbarplot_method_parther_relevance <- plot_grid(
  plot_grid(
    plot.method.time + theme(legend.position = "none"), 
    plot.partner.time + theme(legend.position = "none"),
    plot.conservation.time + theme(legend.position = "none"), 
    ncol = 1, align = "hv", labels = "AUTO")
  , plot_grid(
    get_legend(plot.method.time), 
    get_legend(plot.partner.time),
    get_legend(plot.conservation.time), 
    ncol = 1, align = "hv"), 
  rel_widths = c(1,1,1)
  )
fig_stackedbarplot_method_parther_relevance
```

```{r}
# ggsave(filename = "fig_stackedbarplot_method_partner_relevance.pdf",
       # plot = fig_stackedbarplot_method_parther_relevance,
       # width = 10, height = 10)
```


# Study focus

Make a bar chart of the study focus

```{r bar chart study focus}
lit %>% group_by(Study.Focus) %>% distinct(TI) %>% tally() %>%
 ggplot() +
  geom_bar(aes(reorder(Study.Focus,n),n), stat="identity", width=1, color="white") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("")+
  ylab("Number of Studies")
```


Make a pie chart of the study focus

```{r pie chart study focus}
lit %>% group_by(Study.Focus) %>% distinct(TI) %>% tally() %>%
 ggplot(aes(x="", y=n, fill=Study.Focus)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()
```

# Experimental factors

Make a bar chart of the Experimental Factors

```{r bar chart expt factors}
lit %>% group_by(Expt.Factors) %>% filter(!is.na(Expt.Factors)) %>% distinct(TI) %>% tally() %>%
 ggplot() +
  geom_bar(aes(reorder(Expt.Factors,n),n), stat="identity", width=1, color="white") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("")+
  ylab("Number of Studies")
```

```{r bar chart expt factors filtered}
lit %>% group_by(Expt.Factors) %>% filter(!is.na(Expt.Factors)) %>% distinct(TI) %>% tally() %>%
  filter(n>5) %>%
 ggplot() +
  geom_bar(aes(reorder(Expt.Factors,n),n), stat="identity", width=1, color="white") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("")+
  ylab("Number of Studies")
```


# Focal families and species

Make a bar chart of the Focal Family

```{r bar chart focal family}
lit %>% group_by(Focal.Fam) %>% filter(!is.na(Focal.Fam)) %>% distinct(TI) %>% tally() %>%
 ggplot() +
  geom_bar(aes(reorder(Focal.Fam,n),n), stat="identity", width=1, color="white") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("")+
  ylab("Number of Studies")
```

```{r bar chart focal family filtered}
lit %>% group_by(Focal.Fam) %>% filter(!is.na(Focal.Fam)) %>% distinct(TI) %>% tally() %>%
  filter(n>5) %>%
  ggplot() +
  geom_bar(aes(x=reorder(Focal.Fam,n),y=n), stat="identity", width=1, color="white") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("")+
  ylab("Number of Studies")
```







# Sunburst plots







# Figure 3A: The proportion of experiments by holobiont member, study focus, and experimental factors.

```{r}
lit.sun.2 <- lit %>% select(Partner, Study.Focus, Expt.Factors, TI) %>% distinct()
```

How many experiments? How many publications?

```{r}
nrow(lit.sun.2)
lit.sun.2 %>% group_by(TI) %>% distinct(TI) %>% nrow()
```


Make tallys and percentage by partner

```{r}
lit.sun.2.partner <- lit.sun.2 %>% filter(!is.na(Partner)) %>% 
  group_by(Partner) %>% 
  tally() %>% arrange(desc(n)) %>%
  filter(n>7) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  rename(ids = "Partner") %>%
  mutate(parents = "NA") %>%
  mutate(labels = paste0(ids, "<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```

Make tallys and percentage by study focus for each partner

```{r}
lit.sun.2.study <- lit.sun.2 %>% filter(!is.na(Partner)) %>%
  filter(Partner %in% lit.sun.2.partner$ids) %>%
  group_by(Partner, Study.Focus) %>% 
  tally() %>% arrange(desc(n)) %>%
  # filter(n>4) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  rename(parents = "Partner") %>%
  mutate(ids = paste(Study.Focus,parents,sep="_")) %>%
  mutate(labels = paste0(Study.Focus,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```

Make tallys and percentage by expt factor for each study focus for each partner

```{r}
lit.sun.2.expt <- lit.sun.2 %>% filter(!is.na(Partner)) %>%
  group_by(Partner, Study.Focus, Expt.Factors) %>%
  filter(Partner %in% lit.sun.2.partner$ids) %>%
  tally() %>% arrange(desc(n)) %>%
  # filter(n>0) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  mutate(parents = paste(Study.Focus,Partner,sep="_")) %>%
  mutate(ids = paste(Expt.Factors,parents,sep="_")) %>%
  mutate(labels = paste0(Expt.Factors,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```

```{r}
sunplot.partner2study2expt.data <- bind_rows(lit.sun.2.partner,
                              lit.sun.2.study,
                              lit.sun.2.expt)

colors <- ifelse(grepl("^Coral host$", sunplot.partner2study2expt.data$ids),pal_npg()(1)[1],
                 ifelse(grepl("^Coral host;Sym", sunplot.partner2study2expt.data$ids),pal_npg()(4)[4],
                 ifelse(grepl("^Anemone host$", sunplot.partner2study2expt.data$ids),pal_npg()(5)[5],
                ifelse(grepl("^Anemone host;Sym", sunplot.partner2study2expt.data$ids),pal_npg()(2)[2],
                 ifelse(grepl("^Symbiodiniaceae", sunplot.partner2study2expt.data$ids),pal_npg()(3)[3],
                ifelse(grepl("_Symbiodiniaceae", sunplot.partner2study2expt.data$ids),pal_npg()(7)[7],
                ifelse(grepl("_Coral host$", sunplot.partner2study2expt.data$ids),pal_npg()(1)[1],
                ifelse(grepl("_Coral host;Sym", sunplot.partner2study2expt.data$ids),pal_npg()(4)[4],
                ifelse(grepl("_Anemone host$", sunplot.partner2study2expt.data$ids),pal_npg()(5)[5],
                ifelse(grepl("_Anemone host;Sym", sunplot.partner2study2expt.data$ids),pal_npg()(2)[2],
                ifelse(grepl("^Bacteria", sunplot.partner2study2expt.data$ids),pal_npg()(7)[7],
                ifelse(grepl("_Bacteria", sunplot.partner2study2expt.data$ids),pal_npg()(7)[7],
                 "grey"))))))))))))
sunplot.partner2study2expt.data$colors <- colors
```


Sunplot from study focus to expt factors

```{r}
sunplot.partner2study2expt <- plot_ly(sunplot.partner2study2expt.data, ids = ~ids, labels = ~labels, 
                 parents = ~parents, values= ~n, 
                 type = 'sunburst', branchvalues = 'total',
                 insidetextorientation="auto",
                 width = 500, height = 400,
                 marker = list(colors = sunplot.partner2study2expt.data$colors),
                 color = I('#FFFFFF'))
sunplot.partner2study2expt
```

```{r}
# orca(sunplot.partner2study2expt, file="sunplot_partner2study2expt.pdf")
# orca(sunplot.partner2study2expt, file="sunplot_partner2study2expt.svg")
```









# Figure 3B: The proportion of experiments in coral hosts by life-stage, family, and species.

Make tallys and percentage by developmental stage

```{r}
lit.sun.1 <- lit %>% filter(grepl("host",Partner)) %>%
  select(DevelopmentalStage, FocalSpp, family, TI) %>% distinct() %>% 
  filter(!is.na(DevelopmentalStage)) %>% filter(!DevelopmentalStage=="NA")
```

How many experiments? How many publications?

```{r}
nrow(lit.sun.1)
lit.sun.1 %>% group_by(TI) %>% distinct(TI) %>% nrow()
```

```{r}
lit.sun.1.dev.stage <- lit.sun.1 %>%
  group_by(DevelopmentalStage) %>%
  tally() %>% arrange(desc(n)) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  rename(ids = "DevelopmentalStage") %>%
  mutate(parents = "NA") %>%
  mutate(labels = paste0(ids, "<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```

Make tallys and percentage by focal family for each developmental stage

```{r}
lit.sun.1.focal.fam <- lit.sun.1 %>% 
  filter(DevelopmentalStage %in% lit.sun.1.dev.stage$ids) %>%
  filter(!is.na(family)) %>%
  group_by(DevelopmentalStage, family) %>% 
  tally() %>% arrange(desc(n)) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  rename(parents = "DevelopmentalStage") %>%
  mutate(ids = paste(family,parents,sep="_")) %>%
  mutate(labels = paste0(family,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```


```{r}
sunplot.dev2fam.data <- bind_rows(lit.sun.1.dev.stage,
                              lit.sun.1.focal.fam)

colors <- ifelse(grepl("^adult", sunplot.dev2fam.data$ids),pal_npg()(1)[1],
                 ifelse(grepl("^larvae", sunplot.dev2fam.data$ids),pal_npg()(4)[4],
                 ifelse(grepl("^juvenile", sunplot.dev2fam.data$ids),pal_npg()(2)[2],
                 ifelse(grepl("_adult", sunplot.dev2fam.data$ids),pal_npg()(5)[5],
                 ifelse(grepl("_larvae", sunplot.dev2fam.data$ids),pal_npg()(6)[6],
                 ifelse(grepl("_juvenile", sunplot.dev2fam.data$ids),pal_npg()(2)[2],
                 "grey"))))))

sunplot.dev2fam.data$colors <- colors
sunplot.dev2fam.data$colors <- if_else(grepl("Sym",sunplot.dev2fam.data$ids),
                                       true = pal_npg()(3)[3],
                                       false = sunplot.dev2fam.data$colors)
```

Sunplot from developmental stage to focal family

```{r}
sunplot.dev2fam <- plot_ly(sunplot.dev2fam.data, ids = ~ids, labels = ~labels, 
                 parents = ~parents, values= ~n, 
                 type = 'sunburst', branchvalues = 'total',
                 insidetextorientation="auto",
                 marker = list(colors = sunplot.dev2fam.data$colors),
                 color = I('#FFFFFF'))
sunplot.dev2fam
```

```{r}
# orca(sunplot.dev2fam, file="sunplot_dev2fam.pdf")
# orca(sunplot.dev2fam, file="sunplot_dev2fam.svg")
```


# Developmental Stage --> Focal Family --> Focal Species

Make tallys and percentage by focal species for each focal family for each developmental stage

```{r}
families <- sub("\\_.*", "", lit.sun.1.focal.fam$ids)

lit.sun.1.focal.spp <- lit.sun.1 %>% 
  filter(DevelopmentalStage %in% lit.sun.1.dev.stage$ids) %>% 
  filter(family %in% families) %>%
  group_by(DevelopmentalStage, family, FocalSpp) %>% 
  tally() %>% arrange(desc(n)) %>%
  # filter(n>5) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  mutate(parents = paste(family,DevelopmentalStage,sep="_")) %>%
  mutate(ids = paste(FocalSpp,parents,DevelopmentalStage,sep="_")) %>%
  mutate(labels = paste0(FocalSpp,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```

```{r}
sunplot.dev2fam2spp.data <- bind_rows(lit.sun.1.dev.stage,
                              lit.sun.1.focal.fam,
                              lit.sun.1.focal.spp)

colors <- ifelse(grepl("^adult", sunplot.dev2fam2spp.data$ids),pal_npg()(1)[1],
                 ifelse(grepl("^larvae", sunplot.dev2fam2spp.data$ids),pal_npg()(4)[4],
                 ifelse(grepl("^juvenile", sunplot.dev2fam2spp.data$ids),pal_npg()(2)[2],
                 ifelse(grepl("_adult", sunplot.dev2fam2spp.data$ids),pal_npg()(5)[5],
                 ifelse(grepl("_larvae", sunplot.dev2fam2spp.data$ids),pal_npg()(6)[6],
                 ifelse(grepl("_juvenile", sunplot.dev2fam2spp.data$ids),pal_npg()(2)[2],
                 "grey"))))))

sunplot.dev2fam2spp.data$colors <- colors
sunplot.dev2fam2spp.data$colors <- if_else(grepl("Sym",sunplot.dev2fam2spp.data$ids),
                                       true = pal_npg()(3)[3],
                                       false = sunplot.dev2fam2spp.data$colors)
```

Total number in this chart
```{r}
sunplot.dev2fam2spp.data %>% filter(parents=="NA") %>% summarize(sum(n))
```


Sunplot from developmental stage to focal family to species

```{r}
sunplot.dev2fam2spp <- plot_ly(sunplot.dev2fam2spp.data, ids = ~ids, labels = ~labels, 
                 parents = ~parents, values= ~n, 
                 type = 'sunburst', branchvalues = 'total',
                 insidetextorientation="auto",
                 width = 500, height = 400,
                 marker = list(colors = sunplot.dev2fam2spp.data$colors),
                 color = I('#FFFFFF'))
sunplot.dev2fam2spp
```


```{r}
# orca(sunplot.dev2fam2spp, file="sunplot_dev2fam2spp2.pdf")
# orca(sunplot.dev2fam2spp, file="sunplot_dev2fam2spp2.svg")
```





# Figure 4A: Multistressors

```{r}
lit.sun.6 <- lit %>%
  select(Expt.Factors.Stress, NumStressInt, TI) %>% distinct() %>%
  filter(!is.na(Expt.Factors.Stress)) %>%
  group_by(TI)
```

How many experiments? How many publications?

```{r}
nrow(lit.sun.6)
lit.sun.6 %>% group_by(TI) %>% distinct(TI) %>% nrow()
```

```{r}
lit.sun.6.num <- lit.sun.6 %>% filter(!is.na(NumStressInt)) %>%
  group_by(NumStressInt) %>% 
  tally() %>% arrange(desc(n)) %>%
  filter(n>0) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  mutate(parents = "NA") %>%
  mutate(ids = paste(NumStressInt)) %>%
  mutate(labels = paste0(NumStressInt,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```

```{r}
lit.sun.6.expt.fact<- lit.sun.6 %>%
  filter(!is.na(Expt.Factors.Stress)) %>%
  ungroup() %>%
  group_by(NumStressInt,Expt.Factors.Stress) %>% 
  tally() %>% arrange(desc(n)) %>%
  filter(n>0) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  mutate(parents = paste(NumStressInt)) %>%
  mutate(ids = paste(Expt.Factors.Stress,parents,sep="_")) %>%
  mutate(labels = paste0(Expt.Factors.Stress,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```


```{r}
sunplot.expt2num.data <- bind_rows(lit.sun.6.expt.fact,
                                lit.sun.6.num)

colors <- ifelse(grepl("1", sunplot.expt2num.data$ids),pal_npg()(1)[1],
                 ifelse(grepl("2", sunplot.expt2num.data$ids),pal_npg()(4)[4],
                 ifelse(grepl("3", sunplot.expt2num.data$ids),pal_npg()(3)[3],
                 ifelse(grepl("4", sunplot.expt2num.data$ids),pal_npg()(2)[2],
                 "grey"))))

sunplot.expt2num.data$colors <- colors
```

Sunplot from experimental factors to number of experimental factors

```{r}
sunplot.num2exp <- plot_ly(sunplot.expt2num.data, ids = ~ids, labels = ~labels, 
                 parents = ~parents, values= ~n, 
                 type = 'sunburst', branchvalues = 'total',
                 insidetextorientation="radial",
                 width = 500, height = 400,
                 marker = list(colors = sunplot.expt2num.data$colors),
                 color = I('#FFFFFF'))
sunplot.num2exp
```

```{r}
# orca(sunplot.num2exp, file = "sunplot_num2exp_sub.pdf")
# orca(sunplot.num2exp, file = "sunplot_num2exp_sub.svg")
```







# Figure 3D: The empirical focus of studies with conservation and restoration relevance

```{r}
lit.sun.4 <- lit %>% select(Expt.Factors, Conservation.Relevance, TI) %>% filter(!is.na(Conservation.Relevance)) %>% distinct()
```

How many experiments? How many publications?

```{r}
nrow(lit.sun.4)
lit.sun.4 %>% group_by(TI) %>% distinct(TI) %>% nrow()
94/350
```

Make tallys and percentage by conservation relevance

```{r}
lit.sun.4.cons <- lit.sun.4 %>% filter(!is.na(Conservation.Relevance)) %>%
  group_by(Conservation.Relevance) %>% 
  tally() %>% arrange(desc(n)) %>%
  filter(n>0) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  mutate(parents = "NA") %>%
  mutate(ids = paste(Conservation.Relevance)) %>%
  mutate(labels = paste0(Conservation.Relevance,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```

```{r}
lit.sun.4.expt <- lit.sun.4 %>% filter(!is.na(Expt.Factors)) %>%
  filter(Conservation.Relevance %in% lit.sun.4.cons$ids) %>%
  group_by(Conservation.Relevance,Expt.Factors) %>% 
  tally() %>% arrange(desc(n)) %>%
  filter(n>0) %>%
  mutate(percent = round(n/sum(n),2)*100) %>%
  mutate(parents = paste(Conservation.Relevance)) %>%
  mutate(ids = paste(Expt.Factors,parents,sep="_")) %>%
  mutate(labels = paste0(Expt.Factors,"<br>",percent,"%")) %>%
  select(ids, labels, parents, n, percent)
```


```{r}
sunplot.cons2expt.data <- bind_rows(lit.sun.4.cons,
                                lit.sun.4.expt)

colors <- ifelse(grepl("biomarker", sunplot.cons2expt.data$ids),pal_npg()(1)[1],
                 ifelse(grepl("inter-population variation", sunplot.cons2expt.data$ids),pal_npg()(2)[2],
                 ifelse(grepl("intra-population variation", sunplot.cons2expt.data$ids),pal_npg()(3)[3],
                 ifelse(grepl("management decisions", sunplot.cons2expt.data$ids),pal_npg()(4)[4],
                 ifelse(grepl("reef festoration technology", sunplot.cons2expt.data$ids),pal_npg()(5)[5],
                 ifelse(grepl("heritability", sunplot.cons2expt.data$ids),pal_npg()(6)[6],
                 ifelse(grepl("restoration lines", sunplot.cons2expt.data$ids),pal_npg()(7)[7],
                 ifelse(grepl("cryptic species", sunplot.cons2expt.data$ids),pal_npg()(8)[8],
                 ifelse(grepl("methods development", sunplot.cons2expt.data$ids),pal_npg()(9)[9],
                 "grey")))))))))
sunplot.cons2expt.data$colors <- colors
```

Sunplot from conservation to experimental factors

```{r}
sunplot.cons2expt <- plot_ly(sunplot.cons2expt.data, ids = ~ids, labels = ~labels, 
                 parents = ~parents, values= ~n, 
                 type = 'sunburst', branchvalues = 'total',
                 insidetextorientation="auto",
                 width = 500, height = 400,
                 marker = list(colors = sunplot.cons2expt.data$colors),
                 color = I('#FFFFFF'))
  
sunplot.cons2expt
```

```{r}
# orca(sunplot.cons2expt, file = "sunplot_cons2expt.pdf")
# orca(sunplot.cons2expt, file = "sunplot_cons2expt.svg")
```

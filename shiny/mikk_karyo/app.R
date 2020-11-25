#!/usr/bin/env Rscript

# Set working directory for testing
#library(here)
#setwd(here("shiny", "mikk_karyo"))

# Setup
load_libs = function(){
  library(karyoploteR)
  library(BiocManager)
  options(repos = BiocManager::repositories())
  library(GenomicRanges)
  library(shiny)
  library(dplyr)
}
suppressMessages(load_libs())

# Make custom Karyoplot scaffold

## Get chromosome lengths
med_chr_lens = read.table(file.path("data", "Oryzias_latipes.ASM223467v1.dna.toplevel.fa_chr_counts.txt"),
                          col.names = c("chr", "end"))
## Add start
med_chr_lens$start = 1
## Reorder
med_chr_lens = med_chr_lens %>% 
  dplyr::select(chr, start, end)
## Create custom genome
med_genome = regioneR::toGRanges(med_chr_lens)

# Read in ABBA-BABA sliding windows data and process

## Read in data
df = read.table(file.path("data", "20201022_abba_sliding_windows.txt"), header = T, sep = "\t", as.is = T)
# Convert fd to 0 if D < 0
df$fd = ifelse(df$D < 0,
               0,
               df$fd)
## Change names
df = df %>% 
  dplyr::mutate(p2 = recode(df$p2, hdrr = "HdrR", hni = "HNI", hsok = "HSOK"))
## Set colours
cols <- c("#F3B61F", "#631E68", "#F6673A", "#F33A56", "#55B6B0", "#08605F", "#002642", "#B02156")
names(cols) <- c("HdrR", "HSOK", "HNI", "melastigma", "javanicus", "KW", "HO5", "iCab")
## Filter for melastigma
df_kp = df %>%
  dplyr::filter(p1 == "melastigma")
## make chr numeric
df_kp$scaffold <- as.numeric(df_kp$scaffold)

# Load reference exon density for app

ex_ranges = readRDS(file.path("data", "20201125_exon_ranges.rds"))

# Load MIKK SNP density for app

mikk_ranges = readRDS(file.path("data", "20201125_mikk_ranges.rds"))

# Load exon density for HSOK and HNI

ol_ranges_list = readRDS(file.path("data", "20201125_ol_ranges.rds"))

# Unload packages to save memory

detach("package:BiocManager", unload = T)

# Create app

shiny::shinyApp(
  ui = fluidPage(
    shiny::checkboxGroupInput("chromosome", label = "Chromosome",
                              choices = seq(1, 24), selected = 2),
    shiny::plotOutput("karyoplot")
  ),
  
  server = function(input, output) {
    output$karyoplot = shiny::renderPlot({
      
      kp = plotKaryotype(med_genome, chromosomes = input$chromosome)
      # Add base numbers 
      karyoploteR::kpAddBaseNumbers(kp, tick.dist = 5000000, cex = 0.3)
      # Add data backgrounds
      karyoploteR::kpDataBackground(kp, r0=0, r1 = 1, color = "white")
      # Add axis label
      kpAxis(kp, r0=0.6, r1 = 1, cex = 0.4)
      # Add fd data
      karyoploteR::kpLines(kp,
                           chr = df_kp$scaffold[df_kp$p2 == "HNI"],
                           x = df_kp$mid[df_kp$p2 == "HNI"],
                           y = df_kp$fd[df_kp$p2 == "HNI"],
                           col = "#F6673A",
                           r0=0.6, r1 = 1)
      karyoploteR::kpLines(kp,
                           chr = df_kp$scaffold[df_kp$p2 == "HdrR"],
                           x = df_kp$mid[df_kp$p2 == "HdrR"],
                           y = df_kp$fd[df_kp$p2 == "HdrR"],
                           col = "#F3B61F",
                           r0=0.6, r1 = 1)
      karyoploteR::kpLines(kp,
                           chr = df_kp$scaffold[df_kp$p2 == "HSOK"],
                           x = df_kp$mid[df_kp$p2 == "HSOK"],
                           y = df_kp$fd[df_kp$p2 == "HSOK"],
                           col = "#631E68",
                           r0=0.6, r1 = 1)
      # Add SNP density data
      kpPlotDensity(kp, data=mikk_ranges, col = "#49A379",
                    r0=0, r1=0.2, 
                    window.size = 25000)
      kpPlotDensity(kp, data=ol_ranges_list$hni, col = "#F6673A",
                    r0=0.2, r1=0.4, 
                    window.size = 25000)
      kpPlotDensity(kp, data=ol_ranges_list$hsok, col = "#631E68", 
                    r0=0.4, r1=0.6, 
                    window.size = 25000)
      # Add exon density to ideogram
      kpPlotDensity(kp, data=ex_ranges, col = "#f77cb5",
                    data.panel = "ideogram",
                    window.size = 25000,
                    r0 = 0.5, r1 = 1)
      kpPlotDensity(kp, data=ex_ranges, col = "#f77cb5",
                    data.panel = "ideogram",
                    window.size = 25000,
                    r0 = 0.5, r1 = 0)
      # Add labels
      kpAddLabels(kp, labels="MIKK",
                  r0=0, r1=0.2, 
                  cex = 0.4)
      kpAddLabels(kp, labels="HNI",
                  r0=0.2, r1=0.4, 
                  cex = 0.4)
      kpAddLabels(kp, labels="HSOK",
                  r0=0.4, r1=0.46, 
                  cex = 0.4)
      kpAddLabels(kp, labels=bquote(italic(f[d])),
                  r0=0.6, r1=1, 
                  label.margin = 0.035,
                  cex = 0.6)
      
    })  
  }
)

# Deploy to Shiny.io
#rsconnect::setAccountInfo(name='ian-brettell',
#                          token='33601423D3EB1E5C940DDBD33EABB0CC',
#                          secret='AJtlagvG43M9AsoHl33LaqeQNlO6xXEKaPWqVb5k')
#rsconnect::configureApp("mikk_karyo", size = "large")
#rsconnect::deployApp()



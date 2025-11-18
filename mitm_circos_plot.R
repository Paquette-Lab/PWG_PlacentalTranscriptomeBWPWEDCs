####---------------------------------------------------------------####
## Aim: use PHTHALATE and PLACENTAL EFFICIENCY results to identify 
## meet-in-the-middle (MITM) genes. Then, plot the DEGs and MITM
## genes in a CIRCOS plot using circlize
####---------------------------------------------------------------####

## Load librares
library(tidyverse)
library(here)
library(edgeR)
library(RColorBrewer)
library(circlize)

####---------------------------------------------------------------####
## Load data 

##### full phthalate results
phthalate_all <- readr::read_csv(here::here("Mariana/Draft3/CIRCOS_and_MITMmediation/Phthalate_DEGs_GeoMean_071025.csv")) |>
  dplyr::select(-1)

##### fetal outcomes results 
all_outcome <- readr::read_csv(file = here::here("Mariana", "Env Int Revision", "output", "full_BW_PW_ratio_results.csv"))

## Wrangle data
outcome_results <- all_outcome |> 
    dplyr::mutate(Track_Sublevel = outcome, 
                  SECTOR = outcome,
                  Sector = outcome,
                  Track_Analysis = analysis) |>
  dplyr::select(genes = gene_symbol, 
                logFC = coefficient,
                FDR,
                Sector, 
                SECTOR,
                Track_Sublevel,
                Track_Analysis)

phthalate_results <-   phthalate_all |>
  dplyr::mutate(Track_Analysis = "Overall",
                Sector = "Exposure",
                SECTOR = "Phthalate") |>
  dplyr::select(genes, 
                logFC,
                FDR = adj.P.Val, 
                Sector, 
                SECTOR,
                Track_Sublevel = Phthalate, 
                Metabolite = Phthalate,
                Track_Analysis
  )


circos_data <- dplyr::bind_rows(
  phthalate_results,
  outcome_results
)

####---------------------------------------------------------------####
## Format DATA

Dat <- circos_data |>
  ## remove unneeded columns
  ## note that I tried a couple of layouts, hence the extra Sector and Track
  dplyr::select(-c("Sector", "Track_Sublevel")) |>
  ## filter to only include significant results
  dplyr::filter(FDR<0.05) |>
  dplyr::filter(Track_Analysis == "Overall" & SECTOR != "PW") |>
  ## arrange by sector (chemical class) and metabolite
  dplyr::arrange(SECTOR, Metabolite, Track_Analysis) |> 
  ## convert sector to factor 
  dplyr::mutate(
    Sector = factor(SECTOR, levels = c("BWadj", "BW:PW", "Phthalate"))
  ) |>
  ## assign unique x_position within each sector based on metabolite order
  dplyr::mutate(
    x_position = ave(seq_along(Metabolite), Sector, FUN = seq_along),
    logFC = dplyr::if_else(Sector == "PW", 0, logFC)
  ) 

## first pull out meet in the middle genes
## color in the MEET IN THE MIDDLE results
shared_genes <- intersect(
  Dat$genes[Dat$Sector == "Phthalate"],
  Dat$genes[Dat$Sector %in% c("BWadj", "BW:PW")]
)

####---------------------------------------------------------------####
#### Initialize Plot ####
sector_counts <- table(Dat$Sector)
sector_widths <- as.numeric(sector_counts)

# Initialize circos layout with sector widths
circos.clear()
circos.par(gap.degree = 8)
circos.par(track.margin = c(0, 0))  # Reduce cell padding
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
gap_size <- 1
sector_xlims <- cbind(rep(0, length(sector_counts)), 
                      sector_widths + gap_size)
circos.initialize(factors = names(sector_counts), 
                  xlim = sector_xlims)

####---------------------------------------------------------------####
#### TRACK 1 :: SECTOR #####
sector_colors <- setNames(rep(c("lightsalmon", "lightgoldenrod", "powderblue"),
                              length(levels(Dat$Sector))), levels(Dat$Sector))

# sector_colors <- setNames(rep(c("darkcyan","paleturquoise3","darkseagreen", "seashell"), 
#                               length(levels(Dat$Sector))), levels(Dat$Sector))

# Create the ChemicalClass track
circos.track(ylim = c(0, 0.05), track.height = 0.05, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  if (sector.index %in% names(sector_colors)) {
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 
                0.05, col = sector_colors[sector.index], 
                border = NA)
    # circos.text(mean(CELL_META$xlim), 0.1, sector.index, 
    #             facing = "bending.outside", cex = 1, adj = c(0.5, 1),
    #             niceFacing = TRUE)
  }
})

####---------------------------------------------------------------####
##### TRACK 2 :: METABOLITES #####
# Load color palettes
phthalate_colors <- c("darkseagreen", "darkcyan") # Phthalate colors 
BW_color <- "lightsteelblue"  
# PW_color <- "paleturquoise3"
ratio_color <- "slategray1"  

# Create a lookup table for Metabolite colors based on Sector
metabolite_colors <- setNames(rep("gray", length(unique(Dat$Metabolite))), unique(Dat$Metabolite))

# Assign colors based on the Sector category
for (sector in unique(Dat$Sector)) {
  sector_metabolites <- unique(Dat$Metabolite[Dat$Sector == sector])
  
  if (sector == "Phthalate") {
    metabolite_colors[sector_metabolites] <- phthalate_colors[seq_along(sector_metabolites) %% length(phthalate_colors) + 1]
  } else if (sector == "BW:PW") {
    metabolite_colors[sector_metabolites] <- ratio_color
  } else if (sector == "BWadj") {
    metabolite_colors[sector_metabolites] <- BW_color
  } else if (sector == "PW") {
    metabolite_colors[sector_metabolites] <- PW_color
  }
}

# Create the Metabolite track with labels only on the first row of each metabolite
circos.track(ylim = c(0, 0.05), 
             track.height = 0.05, ## originally 0.05
             panel.fun = function(x, y) {
               sector.index <- CELL_META$sector.index
               metabolites <- Dat[Dat$Sector == sector.index, ]  # Filter metabolites for current sector
               
               if (nrow(metabolites) > 0) {
                 seen_metabolites <- c()  # Keep track of labeled metabolites
                 
                 for (i in seq_len(nrow(metabolites))) {
                   metabolite_name <- metabolites$Metabolite[i]
                   
                   # Draw the colored rectangles
                   circos.rect(metabolites$x_position[i] - 0.5, 0, metabolites$x_position[i] + 0.5, 0.05, 
                               col = metabolite_colors[metabolite_name], 
                               border = NA)
                   
                   #   # Only label the first occurrence of each metabolite
                   #     if (!(metabolite_name %in% seen_metabolites)) {
                   #   circos.text(x = metabolites$x_position[i], y = 0,
                   #   labels = metabolite_name,
                   #   niceFacing = TRUE,
                   #   facing = "clockwise", cex = 0.5, adj = c(0, 1))
                   #   seen_metabolites <- c(seen_metabolites, metabolite_name)  # Mark as labeled
                   #   }
                   
                 }
               }
             })

####---------------------------------------------------------------####
##### TRACK 3 ######

# Function to map logFC values to colors
map_colors <- function(values) {
  pal_pos <- colorRampPalette(c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)  # Light pink to red
  pal_neg <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0"))(100)  # Light blue to blue
  colors <- rep("white", length(values))
  colors[values > 0] <- pal_pos[findInterval(values[values > 0], seq(0, max(values), length.out = 101), all.inside = TRUE)]
  colors[values < 0] <- pal_neg[findInterval(values[values < 0], seq(min(values), 0, length.out = 101), all.inside = TRUE)]
  return(colors)
}

# Apply the color mapping function
Dat$Color <- map_colors(Dat$logFC)

# Identify genes that appear in at least two different sectors
gene_sector_counts <- aggregate(Sector ~ genes, data = Dat, FUN = function(x) length(unique(x)))
shared_across_sectors <- gene_sector_counts$genes[gene_sector_counts$Sector > 1]  # Only keep genes in multiple sectors

# sector_y_limits <- rbind(
#   ## Birthweight 
#   c(min(Dat[Dat$Sector == "BWadj", ]$logFC) * 1.05, max(Dat[Dat$Sector == "BWadj", ]$logFC)*1.05),
#   ## BW: PW
#   c(min(Dat[Dat$Sector == "BW:PW", ]$logFC) * 1.05, max(Dat[Dat$Sector == "BW:PW", ]$logFC)*1.05),
#   ## Phthalate
#   c(min(Dat[Dat$Sector == "Phthalate", ]$logFC) * 1.05, max(Dat[Dat$Sector == "Phthalate", ]$logFC)*1.05)
# )

# Create the heatmap track with only cross-sector shared gene labels
circos.track(sectors = levels(Dat$Sector),
             ylim = c(1.05*min(Dat$logFC), 1.05*max(Dat$logFC)), 
             panel.fun = function(x, y) {
               sector.index <- CELL_META$sector.index
               sector_data <- Dat[Dat$Sector == sector.index, ]
               
               ## axis
               circos.yaxis(side = "right",
                            labels.cex = 0.5,
                            lwd = 0.1,
                            labels.niceFacing = TRUE)
               
               for (i in seq_len(nrow(sector_data))) {
                 gene_name <- sector_data$genes[i]
                 
                 
                 ## Draw heatmap bars
                 circos.rect(sector_data$x_position[i] - 0.5, 0,
                             sector_data$x_position[i] + 0.5, sector_data$logFC[i],
                             col = sector_data$Color[i], border = NA)
                 
                 #  # Label only genes that are **shared across different sectors**
                 if (gene_name %in% shared_genes) {
                   
                   circos.text(sector_data$x_position[i],
                               max(sector_data$logFC),
                               gene_name,
                               facing = "clockwise",
                               niceFacing = TRUE,
                               cex = 0.5,
                               adj = c(1, 0.5)
                   )
                 }
               }
             })

####---------------------------------------------------------------####
####### Interior #######
# Get unique genes in the dataset
unique_genes <- unique(Dat$genes)

# Create a data frame to store link connections
link_data <- data.frame(from_sector = character(), from_x = numeric(),
                        to_sector = character(), to_x = numeric(), gene = character(),
                        stringsAsFactors = FALSE)

# Loop through each unique gene and find its occurrences
for (gene in unique_genes) {
  gene_rows <- Dat[Dat$genes == gene, ] 
  if (nrow(gene_rows) > 1) {
    for (i in 1:(nrow(gene_rows) - 1)) {
      for (j in (i + 1):nrow(gene_rows)) {
        link_data <- rbind(link_data, data.frame(
          from_sector = gene_rows$Sector[i], from_x = gene_rows$x_position[i],
          to_sector = gene_rows$Sector[j], to_x = gene_rows$x_position[j],
          gene = gene
        ))
      }
    }
  } 
}
## meet in the middle
## color in the MEET IN THE MIDDLE results
shared_genes <- intersect(
  Dat$genes[Dat$Sector == "Phthalate"],
  Dat$genes[Dat$Sector %in% c("BWadj", "BW:PW", "PW")]
)

mitm_a <- Dat |>
  dplyr::filter(Sector == "Phthalate") |>
  dplyr::mutate(phthalate_x = x_position) |>
  dplyr::select(genes, Metabolite, phthalate_x) |>
  ## get colors by linking to metabolite colors
  dplyr::left_join(enframe(metabolite_colors) |> 
                     # dplyr::filter(stringr::str_detect(name, "^M")) |>
                     dplyr::select(Metabolite = name, color = value))

mitm_b <- Dat |>
  dplyr::filter(Sector != "Phthalate") |>
  dplyr::mutate(outcome_x = x_position) |>
  dplyr::select(genes, Sector, outcome_x)

mitm <- full_join(mitm_a, mitm_b) |>
  dplyr::filter(genes %in% shared_genes)

# Draw the links between phthalates using circlize
for (i in 1:nrow(link_data)) {
  from_sector <- link_data$from_sector[i]
  to_sector <- link_data$to_sector[i]
  from_x <- link_data$from_x[i]
  to_x <- link_data$to_x[i]
  gene_name <- link_data$gene[i]
  
  metab_row <- Dat[Dat$genes == gene_name & Dat$Sector == "Phthalate", ]
  metab <- if (nrow(metab_row) > 0) metab_row$Metabolite[1] else NA
  link_col <- NA
  link_lwd <- 1
  
  if (from_sector == "Phthalate" && to_sector == "Phthalate") {
    link_col <- "lightgrey"
    link_lwd <- 1
  } else if (from_sector == "Phthalate" && to_sector %in% c("BWadj", "BW:PW", "PW")) {
    if (!is.na(metab) && metab %in% names(metabolite_colors)) {
      if (metab %in% names(metabolite_colors)) {
        link_col <- metabolite_colors[metab]
      } else {
        cat("  --> Metabolite color missing for:", metab, "\n")
        link_col <- "orange"
      }
      link_lwd <- 2
    }
  } else {
    next
  }
  
  if (!is.na(link_col)) {
    circos.link(from_sector, from_x, to_sector, to_x, col = link_col, border = NA, lwd = link_lwd)
  }
}

## color in the MEET IN THE MIDDLE results
for (i in 1:nrow(mitm)) {
  
  circos.link("Phthalate", mitm$phthalate_x[i],
              mitm$Sector[i], mitm$outcome_x[i],
              col = mitm$color[i], lwd = 2)
}

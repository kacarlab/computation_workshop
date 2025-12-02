# growthcurver_analysis.R
# Analyze growth curve data using R package Growthcurver
# See https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html for more information on Growthcurver package
# Updated: 2023-06-12, Amanda K. Garcia (akgarcia3@wisc.edu)


### CHANGE THESE VARIABLES AS NEEDED ###

# CSV growth data file path
data_fname <- "/Users/fer/Desktop/week4/IF2_NTD_42C_cleaned.csv"

# Specify output location
out_fpath <- "/Users/fer/Desktop/week4/"

# Replace with "hours, seconds", or "minutes" (check spelling). 
# Will automatically convert to hours
time_units <- "minutes"

# Specify cutoff time point *in hours* to trim data.
# Cutoff is inclusive (e.g., "24" will include first 24 hours)
# Value of "0" means no trimming
# Value that is outside of data time range will result in no trimming
trim_at_time <-0

# Define sample groups (e.g., strains or conditions)
sample_group_names <- c("WT-IF2","dN-IF2", "dG1-IF2")

# Define columns in CSV growth data file that belong to each group
# Numbering starts AFTER the "time" column
# E.g., column 1 is first sample column with measurement data
# List order must match sample group names above
sample_group_cols <- list(1:5,6:10,11:15)

# NOTE: Will automatically perform background correction if a "blank" column is provided.
# Values in blank column are subtracted from all other columns per row (bg_correct = "blank")
# If "blank" column is not provided, will correct background by subtracting all values in a column
# by the minimum value in the same column (bg_correct = "min").
# This method is suitable if the background is NOT expected to change during the experiment.

#############


### SETUP ###
#install.packages(c("growthcurver","dplyr","gridExtra","ggplot2","pals", "tidyr","ggpubr","pdftools","grid","tools"))
#Load libraries
library(growthcurver)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(pals)
library(tidyr)
library(ggpubr)
library(pdftools)
library(grid)
library(tools)


### IMPORT GROWTH DATA ###

# Import data
data <- read.csv(file = data_fname, header = TRUE, sep = ",")


# Define output file basename
out_base <- paste(out_fpath, file_path_sans_ext(basename(data_fname)), sep="")


### CLEAN UP DATA ###

# Convert time units if provided in minutes or seconds
if (time_units=="minutes") {
  print("Converting time units from minutes to hours")
  data$time <- data$time / 60 # Convert minutes to hours
} else if (time_units=="seconds") {
  print("Converting time units from seconds to hours")
  data$time <- data$time / 60 / 60 # Convert seconds to hours
} else if (time_units=="hours") {
  print("Time units are in hours")
} else {
  print("Cannot recognize specified time units. Will assume in hours")
}

# Clean up na values
data[is.na(data)] <- 0 

# Trim data (optional)

if (trim_at_time != 0 & trim_at_time <= (tail(data$time, n=1)-1)) {
  data <- data %>% filter(row_number() < which(time==(trim_at_time+1)))
  print(paste("Data trimmed to", trim_at_time, "hours"))
} else {
  print("Warning: Data not trimmed because trim_at_time value is either 0 or outside data time range")
}

### SUMMARIZE DATA ###

# Summarize data and export fitted growth curves

if ("blank" %in% colnames(data)) {
  data_summary_prelim <- SummarizeGrowthByPlate(
    plate=data,
    bg_correct = "blank",
    plot_fit = TRUE,
    plot_file = paste(out_base, "_fitted_curves.pdf", sep=""))
  data_summary <- head(data_summary_prelim, -1) # Remove residual blank row from data summary
} else {
  data_summary <- SummarizeGrowthByPlate(
    plate=data,
    bg_correct = "min",
    plot_fit = TRUE,
    plot_file = paste(out_base, "_fitted_curves.pdf", sep=""))
} 

# IMPORTANT to inspect fitted growth curve PDF for samples with poor fit
# Poor fit might result from poorly resolved lag and stationary phase
# May be necessary to trim data to improve fit

### CALCULATE SAMPLE GROUP STATISTICS ###

# Add group labels to data summary
data_summary$group <- NA

for (i in 1:(length(sample_group_names))) {
  data_summary$group[unlist(sample_group_cols[i])] <- sample_group_names[i]
}  

# Calculate summary statistics by sample group
data_summary_stats <- data_summary %>%
  group_by(group) %>%
  summarise(
    mean_t_gen=mean(t_gen),
    sd_t_gen=sd(t_gen),
    mean_t_mid=mean(t_mid),
    sd_t_mid=sd(t_mid)
  )
data_summary_stats <- data_summary_stats[match(sample_group_names, data_summary_stats$group),] # Change group ordering


### PLOT DATA ###

# Make barplots of growth parameters
data_summary_stats$group <- factor(data_summary_stats$group,  # Change group ordering
                                   levels = sample_group_names)

# Doubling time plot (t_gen)
t_gen_barplot <- ggplot(data=data_summary_stats, aes(x=group, level=sample_group_names, y=mean_t_gen, fill=group)) +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_t_gen-sd_t_gen, ymax=mean_t_gen+sd_t_gen), width=.05) +
  labs(x="Group", y="Doubling time (h)") +
  scale_fill_manual(name = "Group", values = cols25()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Midpoint time plot (t_mid)
t_mid_barplot <- ggplot(data=data_summary_stats, aes(x=group, y=mean_t_mid, fill=group)) +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_t_mid-sd_t_mid, ymax=mean_t_mid+sd_t_mid), width=.05) +
  labs(x="Group", y="Midpoint time (h)") +
  scale_fill_manual(name = "Group", values = cols25()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Make growth curve scatter plot

# If "blank" provided, subtract blank value from all measurements
if ("blank" %in% colnames(data)) {
  data_blanked <- data
  data_blanked[2:ncol(data_blanked)] <- data_blanked[2:ncol(data_blanked)]-data_blanked$blank
  data_blanked <- data_blanked[,!names(data_blanked) %in% c("blank")]
} else {
  data_blanked <- data # do nothing
}

# Reformat data for scatter plot
## so I think averaging the 6 replicates of each biological rep to make
## total of 3 data points (1 per bio rep) per time point per set is good.
data_for_scatterplt <- data_blanked %>% gather(sample, OD, -time)

data_for_scatterplt$group <- NA

for (i in 1:(length(sample_group_names))) {
  rows = c()
  for (j in sample_group_cols[[i]]) {
    new_rows = (j*nrow(data)-(nrow(data)-1)):(j*nrow(data))
    rows = append(rows, new_rows)
  }
  data_for_scatterplt$group[rows] <- sample_group_names[i]
} 

# Make scatter plot
y_ticks <- -10:10
y_ticks_log <- sapply(y_ticks, function(x) exp(x))

growth_curve_plot <- ggplot(data_for_scatterplt, aes(x = time, y = OD)) +
  geom_point(alpha = 0.35, aes(color = group), size = 0.5) +
  scale_y_continuous(trans = "log", breaks = y_ticks_log, labels = y_ticks) +
  scale_color_manual(name = "Group", values = cols25(), breaks = sample_group_names) +
  labs(x="Time (h)", y="ln OD") +
  theme_classic()

# Add smoothed line per sample group
for (i in 1:(length(sample_group_names))) {
  growth_curve_plot <- growth_curve_plot + geom_smooth(
    data = data_for_scatterplt[data_for_scatterplt$group==(sample_group_names[i]),],
    se = TRUE, color = (cols25())[i])
}


### SAVE REPORT ###

# Open PDF device
pdf(file=(paste(out_fpath, "temp_plots.pdf", sep="")), height=8.5, width=11)

# Save summary table
total_rows_per_page = 38 
start_row = 1 

if(total_rows_per_page > nrow(data_summary)){
  end_row = nrow(data_summary)
}else {
  end_row = total_rows_per_page 
}    

for(i in 1:ceiling(nrow(data_summary)/total_rows_per_page)){
  
  grid.newpage()   
  
  grid.table(data_summary[start_row:end_row, ],
             theme=ttheme_minimal(
               core=list(fg_params=list(cex=0.5)),
               colhead=list(fg_params=list(cex=0.5)),
               rowhead=list(fg_params=list(cex=0.5)))
  )             
  
  start_row = end_row + 1
  
  if((total_rows_per_page + end_row) < nrow(data_summary)){
    
    end_row = total_rows_per_page + end_row
    
  }else {
    
    end_row = nrow(data_summary)
  }    
}

# Save summary stats table
grid.newpage()
grid.table(data_summary_stats, theme=ttheme_minimal(
  core=list(fg_params=list(cex=0.5)),
  colhead=list(fg_params=list(cex=0.5)),
  rowhead=list(fg_params=list(cex=0.5))
))

# Save plots
ggarrange(growth_curve_plot, 
          ggarrange(t_gen_barplot, t_mid_barplot, labels = c("B", "C"), ncol = 2),
          labels = "A",
          nrow = 2
)

dev.off()

# Combine with fitted curve PDF into single report
pdf_combine(c(paste(out_base, "_fitted_curves.pdf", sep=""), paste(out_fpath, "temp_plots.pdf", sep="")),
            output=(paste(out_base, "_report.pdf", sep="")))
file.remove(c(paste(out_base, "_fitted_curves.pdf", sep=""), paste(out_fpath, "temp_plots.pdf", sep="")))

# Write summary and summary stats to csv files
write.csv(data_summary, file = (paste(out_base, "_summary.csv", sep="")), row.names=FALSE)
write.csv(data_summary_stats, file = (paste(out_base, "_summary_stats.csv", sep="")), row.names = FALSE, quote = FALSE)

print(paste("Saved report to ", out_base, "_report.pdf", sep=""))
print(paste("Saved data summary files to ", out_base, "_summary.csv and ",  out_base, "_summary_stats.csv", sep=""))

print("Finished with Growthcurver analysis")


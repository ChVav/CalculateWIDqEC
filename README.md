# CalculateWIDqEC v0.0.1

The purpose of the WID-qEC is to estimate the risk for endometrial (WID-qEC) or cervical cancer (WID-qCIN) based on the methylation status of DNA recovered from cervical smear samples, assayed with MethyLight.

The CalculateWIDqEC App needs the exported "Results file" from Quantstudio, and assumes the export has been done using the correct sampleIDs, targets and fixed Ct thresholds.

## Instructions

Make sure R is installed and in your environment path.

Make sure your default browser is a modern one.

Libraries needed: shiny, dplyr, tidyverse, openxlsx, ggplot2, ggpubr.

Double-click CalculateWIDqEC.bat (or a short-cut thereof), and follow instructions in the browser.
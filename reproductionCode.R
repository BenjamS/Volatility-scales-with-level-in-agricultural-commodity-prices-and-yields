# ==========================================================================
# Title: Reproduction code for "Fluctuation scaling in agricultural commodity markets"
# Author: Benjamin Schiek (Alliance of Bioversity International and CIAT)
# Date: April 2026
# Journal: PNAS Nexus (Brief Report)
# Description: This script imports raw FAO and World Bank data, performs deflation (MUV Index),
#              first-differencing, fits the models ln(sd(price or yield))~ln(mean(price or yield)),
#              performs Shapiro Wilk test for (log-)normality, and generates 
#              the scaling plots appearing in Figure 1 and the dispersion stability plots
#              in Figure 2.
#===========================================================================
library(tidyverse)
library(patchwork)
# Plot parameters
shapeVec <- c(21:25, 4)
#shapeVec <- c(20:25, 4)
pointSize <- 1.5
smallPointSize <- 1
labelSize <- 2.5
smallLabelSize <- 2
titleSize <- 8
subtitleSize <- 7
legendTextSize <- 6
axisTextSize <- 6
axisTitleSize <- 7
facetTitleSize <- 7
#Folder
thisFolder <- "...your_folder/where_you_put/the_raw_data_files/"
thisFolder <- "D:/OneDrive - CGIAR/Documents 1/CIAT 2/Volatility is proportionate to level in agricultural commodity prices/Data/"
#===========================================================================
# Define utility functions
processDf <- function(thisFilepath, theseElements, theseUnits, isExpPrice = F){
  dfRaw <- read.csv(thisFilepath, stringsAsFactors = F)
  dfRaw <- dfRaw %>% subset(Element %in% theseElements &
                              Area == "World" & Unit %in% theseUnits) %>%
    select(-contains(".")) %>% select(-matches("[NF]$")) %>% rename_with(~ str_remove(., "^Y"))
  gathercols <- colnames(dfRaw)[-c(1:4)]
  df <- dfRaw %>%
    select(-Unit) %>%
    pivot_longer(cols = all_of(gathercols), names_to = "year", values_to = "val")
  # Keep only items with complete time series
  df <- df %>% group_by(Area, Element, Item) %>% mutate(nYrs = length(year))
  #hist(df$nYrs)
  df <- df %>% subset(nYrs == 64)
  # Keep only items with no zero values in their series
  df <- df %>% group_by(Element, Item) %>% mutate(n0 = sum(val == 0)) %>% subset(n0 == 0)
  # Keep only items with no NA values in their series
  df <- df %>% group_by(Element, Item) %>% mutate(nNA = sum(is.na(val))) %>% subset(nNA == 0)
  if(isExpPrice){
    df <- df %>% pivot_wider(names_from = Element, values_from = val) %>%
      mutate(price = `Export value` / `Export quantity`) %>%
      filter(is.finite(price))
  }
  return(df)
}
#----------------------------------------------------------------------------
categorizeItems <- function(thisDf){
  outDf <- thisDf %>%
    mutate(group = case_when(
      # 1. Animal products
      str_detect(Item, regex("Animal|Dairy|Meat|Bovine|Beef|Pig|Poultry|Sheep|Goat|Milk|Cheese|Butter|Egg|Yoghurt|Whey|Fat of|Honey|Beeswax|Tallow|Ghee|Offal|Lard|Fat|Fatty|Cream|Yoghurt", ignore_case = TRUE)) ~ "Animal products",
      # 2. Cereals & starchy roots
      str_detect(Item, regex("Chicory|Taro|Wheat|Fonio|Rice|Maize|Barley|Oats|Sorghum|Millet|Cereal|Buckwheat|Rye|Grain|Bran|Cassava|Potato|Sweet potatoes|Yams|Roots and Tubers|Starch of|Bread|Pastry|pasta", ignore_case = TRUE)) ~ "Cereals & starchy roots",
      # 3. Fruit & vegetables
      str_detect(Item, regex("Beet pulp|Mangoes|Cucumber|Eggplants (aubergines)|Pumpkins|Lemons|Okra|Papayas|Lettuce|Vegetable|Cherries|Apples|Bananas|Oranges|Grapes|Citrus|Berries|Onions|Tomatoes|Vegetables|Fruit|Juice|Melons|Mushrooms|Asparagus|Cabbages|Pineapples|Plantains|Garlic|Spinach|Peaches|Pears|Apricots|Raisins|Artichoke|Avocado|Cauliflower|Chillies|Carrot|Currants|Dates|Figs|Cider|Plums|Persimmons|Tangerines|Quinces", ignore_case = TRUE)) ~ "Fruit & vegetables",
      # 4. Pulses, nuts, oils
      str_detect(Item, regex("Beans|palm fruit|Peas|Lentils|Chick peas|Pulses|Almonds|Cashew|Walnut|Hazelnut|Nuts|Pistachio|Areca|Oil|Soy|Sunflower|Rape|Palm|Coconut|Olive|Linseed|Sesame|Cake of|Groundnut|Mustard|Castor|Margarine|Copra|Karite|seed", ignore_case = TRUE)) ~ "Pulses, nuts, oils",
      # 5. Sugar & stimulants
      str_detect(Item, regex("Anise|Spice|Hop|Sugar|stimulant|Molasses|Coffee|Tea|Cocoa|Maté|Confectionery|Chocolate|Glucose|Lactose|Alcoholic|Wine|Beer|Beverages|Ethyl alcohol|Cinnamon|Cloves|Ginger|Nutmeg|Malt|Pepper|Vanilla", ignore_case = TRUE)) ~ "Sugar, spices, stimulants",
      # 6. Industrial/non-food
      str_detect(Item, regex("Pyrethrum, dried flowers|Ramie, raw or retted|Sisal, raw|Pellets|Cotton|Wool|Silk|Rubber|Tobacco|Hides|Skins|Abaca|Jute|Coir|Flax|Cigarettes|Cigars|Straw|Hair|Waste|Residues|wax|feed|Hay|forage|materials|Textile|Fibre|cat food", ignore_case = TRUE)) ~ "Industrial/feed/non-food",
      # Default for miscellaneous prep or missed items
      TRUE ~ "Other"
    ))
  return(outDf)
}
#----------------------------------------------------------------------------
processRiskRewardDf <- function(thisDf){
  dfRiskReward <- thisDf %>% group_by(group, Item) %>% summarise(mlPrice = mean(log(priceSeries)),
                                                                 slPrice = sd(log(priceSeries)),
                                                                 mPrice = mean(priceSeries),
                                                                 sPrice = sd(diff(priceSeries)[-1]),
                                                                 mDifLprice = mean(diff(log(priceSeries))[-1]),
                                                                 sDifLprice = sd(diff(log(priceSeries))[-1])) %>%
    mutate(lmPrice = log(mPrice), lsPrice = log(sPrice)) %>% select(Item, lmPrice, lsPrice, mlPrice, slPrice,
                                                                    mDifLprice, sDifLprice, group)
  return(dfRiskReward)
}
#----------------------------------------------------------------------------
# lm() wrapper to get the slope, y-intercept, and adj. R-squared reported in
# the top of the panels in Figure 1
getSlopeEtc <- function(dfRiskReward){
  dfMod <- dfRiskReward %>% as.data.frame() %>% select(lsPrice, lmPrice)
  row.names(dfMod) <- dfRiskReward$Item
  mod <- lm(lsPrice ~., dfMod)
  #summary(mod)
  b <- round(mod$coefficients[1], 2)
  m <- round(mod$coefficients[2], 2)
  dfOut <- as.data.frame(broom::glance(mod))
  adjR2 <- round(dfOut$adj.r.squared, 2)
  N <- nobs(mod)
  outVec <- c(slope = m, yInt = b, adjR2 = adjR2, N = N)
  return(outVec)
}
#============================================================================
# End utility function definitions
#============================================================================
#============================================================================
#============================================================================
# Figure 1, panel A - Fluctuation scaling in export price series
#============================================================================
# Import and process the raw FAO export value and quantity data
# and World Bank MUV series for deflation of nominal dollar values.
# Notes: FAO uses the MUV series to deflate Value of Agricultural Production (VAP)
# Here we are doing the same for export values.
# Raw FAO export data taken from FAOSTAT: https://www.fao.org/faostat/en/#data/TCL ("All Data" bulk download)
# Annual WB MUV taken from the "Pink Sheet": https://www.worldbank.org/en/research/commodity-markets
#----------------------------------------------------------------------------
# Import and basic processing of FAO export data
thisFile <- "Trade_CropsLivestock_E_All_Data_exportOnly.csv"
thisFilepath <- paste0(thisFolder, thisFile)
theseElements <- c("Export quantity", "Export value")
theseUnits <- c("t", "1000 USD")
df <- processDf(thisFilepath, theseElements, theseUnits, isExpPrice = T)
# Import and basic processing of WB MUV index
thisFile <- "CMO-Historical-Data-Annual.xlsx"
thisFilepath <- paste0(thisFolder, thisFile)
dfMUV <- readxl::read_xlsx(thisFilepath, sheet = "Annual Indices (Real)")
dfMUV <- dfMUV %>% select(year = 1, MUV = 18) %>% slice(-c(1:8))
dfMUV$MUV <- as.numeric(dfMUV$MUV)
# Rebase MUV to more suitable year dollars (default is 2010 dollars)
# (Rebase to match FAO 2014-2016 constant USD used in FAO VAP data)
muvRebaser <- dfMUV %>%
  filter(year %in% c("2014", "2015", "2016")) %>%
  summarise(muvRebaser = mean(MUV)) %>%
  pull(muvRebaser)
df <- df %>% merge(dfMUV) %>% mutate(priceSeries = price * muvRebaser / MUV)
#----------------------------------------------------------------------------
# Quick Shapiro-Wilk test for lognormality in export price series
# (Mentioned in the Discussion)
dfTest <- df %>% group_by(Item) %>% mutate(shapiro_p = shapiro.test(log(priceSeries))$p.value) %>%
  mutate(stndev = sd(diff(log(priceSeries))[-1])) %>% mutate(m = mean(diff(log(priceSeries))[-1])) %>%
  ungroup() %>% subset(year == 2024) %>% select(Item, shapiro_p, stndev, m)
nShapiroExp <- dfTest %>% subset(shapiro_p >= 0.05) %>% nrow
pctShapExp <- nShapiroExp / nrow(dfTest) * 100
print(paste("Num. export price series passing Shapiro-Wilk:", nShapiroExp))
print(paste("Pct. export price series passing Shapiro-Wilk:", round(pctShapExp, 1)))
#----------------------------------------------------------------------------
# Add commodity group column
df <- categorizeItems(df)
# Take means and standard deviations of the series
dfRiskReward <- processRiskRewardDf(df)
# Drop a few outliers
dfRiskReward <- dfRiskReward %>% subset(group != "Other" & Item != "Crude Materials nes")
# Save for use in Figure 2
dfRRexp <- dfRiskReward; dfRRexp$facetKey <- "A) Export price (constant 2014-16 USD/kg)"
# Set group colors to be used in all Figure 1 plots
n <- length(unique(dfRiskReward$group))
bag_of_colors <- randomcoloR::distinctColorPalette(k = 2 * n)
colorVec <- sample(bag_of_colors, n)
# Get slope and intercept of the fluctuation scaling and add to panel title
outVec <- getSlopeEtc(dfRiskReward)
m <- outVec[1]; b <- outVec[2]; adjR2 <- outVec[3]; N <- outVec[4];
thisFacet <- paste0("N = ", N, " agricultural commodities", "\nAdj. R-squared = ", adjR2, ", Slope = ", m, ", Y intercept = ", b)
dfRiskReward$facetKey <- paste0("A) Export price (constant 2014-16 USD/kg)\n", thisFacet)
#----------------------------------------------------------------------------
# Label items in the export price plot
labelItems <- c("Maize (corn)", "Flour of maize",
                #"Wheat", "Flour, wheat",
                #"Starch, cassava",
                "Coffee, green", "Coffee, decaffeinated or roasted", #"Bananas",
                #"Plantains",
                #"Sweet potatoes",
                #"Potatoes",
                "Cocoa beans",
                #"Wool, greasy", 
                #"Beans, dry",
                "Raw silk (not thrown)",
                #"Apples",
                "Vanilla, raw",
                "Soya beans", "Soya bean oil", "Natural Rubber",# "Rice",
                "Cocoa butter, fat and oil")
#----------------------------------------------------------------------------
# Plot panel A of Figure 1
dfPlot1 <- dfRiskReward
dfPlot1$labelThese <- NA; u <- dfPlot1$Item
dfPlot1$labelThese[which(u %in% labelItems)] <- u[which(u %in% labelItems)]
dfPlot1$group <- as.factor(dfPlot1$group)
names(colorVec) <- levels(dfPlot1$group)
names(shapeVec) <- levels(dfPlot1$group)
fillScale <- scale_fill_manual(name = "group", values = colorVec)
shapeScale <- scale_shape_manual(name = "group", values = shapeVec)
gg <- ggplot(dfPlot1, aes(x = lmPrice, y = lsPrice,
                          group = group,
                          fill = group,
                          shape = group))
gg <- gg + geom_smooth(aes(group = NULL, fill = NULL, shape = NULL), method = lm, se = F)
gg <- gg + geom_point(alpha = 0.6, size = pointSize, color = "black", stroke = 0.5)
gg <- gg + shapeScale
gg <- gg + fillScale
gg <- gg + labs(x = "Nat. log of mean", y = "Nat. log of stand. dev.")
gg <- gg + ggrepel::geom_text_repel(aes(label = labelThese), size = labelSize,
                                    point.padding = 1,
                                    box.padding = 1,
                                    nudge_x = 1,
                                    nudge_y = 1,
                                    max.overlaps = Inf)
gg <- gg + facet_wrap(~facetKey)
gg <- gg + theme_bw()
gg <- gg + theme(axis.title = element_text(size = axisTitleSize),
                 axis.text = element_text(size = axisTextSize),
                 legend.position = "none",
                 strip.background = element_rect(fill = "white"),
                 strip.text = element_text(hjust = 0, size = axisTitleSize))
ggFig1A <- gg
#============================================================================
# Figure 1, panels B-D - fluctuation scaling in producer price and yield series
#============================================================================
# Import and process the raw FAO production and VAP data
# Notes: Production quantity and area harvested data downloaded from
# https://www.fao.org/faostat/en/#data/QCL ("All Data" bulk download)
# Value of Agricultural Production (VAP) data downloaded from
# https://www.fao.org/faostat/en/#data/QV ("All Data" bulk download)
# VAP data already deflated by the World Bank MUV.
# Producer price per kg calculated as deflated VAP / production.
# Producer price per hectare calculated as deflated VAP / area harvested.
thisFile <- "Production_Crops_Livestock_E_All_Data.csv"
thisFilepath <- paste0(thisFolder, thisFile)
theseElements <- c("Production", "Area harvested", "Yield"); theseUnits <- c("t", "ha", "kg/ha")
dfAreaProd <- processDf(thisFilepath, theseElements, theseUnits, isExpPrice = F) %>%
  spread(Element, val) %>% select(Item, year, Production, `Area harvested`, Yield)
thisFile <- "Value_of_Production_E_All_Data.csv"
thisFilepath <- paste0(thisFolder, thisFile)
theseElements <- "Gross Production Value (constant 2014-2016 thousand US$)"; theseUnits <- "1000 USD"
dfVAP <- processDf(thisFilepath, theseElements, theseUnits, isExpPrice = F) %>%
  select(Item, year, val) %>% rename(VAP = val); dfVAP$Element <- NULL
# Calculate kg and hectare producer prices
df <- dfAreaProd %>% merge(dfVAP, by = c("Item", "year")) %>%
  mutate(haPrice = 1000 * VAP / `Area harvested`, kgPrice = VAP / Production,
         year = as.integer(year))
#----------------------------------------------------------------------------
# Add commodity group column
df <- categorizeItems(df)
#----------------------------------------------------------------------------
# Shapiro-Wilk test for (log-)normality in producer price series (mentioned in the Discussion)
dfTest <- df %>% group_by(Item) %>% mutate(nYrs = sum(!is.na(kgPrice))) %>% subset(nYrs == 64) %>%
  mutate(shapiro_p = shapiro.test(scale(log(kgPrice)))$p.value) %>%
  ungroup() %>% subset(year == 2024) %>% select(Item, shapiro_p)
nShapiroKg <- dfTest %>% subset(shapiro_p >= 0.05) %>% nrow
pctShapKg <- nShapiroKg / (dfTest %>% nrow) * 100
print(paste("Num. prod. price/kg series passing Shapiro-Wilk:", nShapiroKg))
print(paste("Pct. prod. price/kg series passing Shapiro-Wilk:", round(pctShapKg, 1)))
dfTest <- df %>% group_by(Item) %>% mutate(nYrs = sum(!is.na(haPrice))) %>% subset(nYrs == 64) %>%
  mutate(shapiro_p = shapiro.test(log(haPrice))$p.value) %>%
  ungroup() %>% subset(year == 2024) %>% select(Item, shapiro_p)
nShapiroHa <- dfTest %>% subset(shapiro_p >= 0.05) %>% nrow
pctShapHa <- nShapiroHa / (dfTest %>% nrow) * 100
print(paste("Num. prod. price/ha series passing Shapiro-Wilk:", nShapiroHa))
print(paste("Pct. prod. price/ha series passing Shapiro-Wilk:", round(pctShapHa, 1)))
#----------------------------------------------------------------------------
# Subset the data for the prodcer price/kg plot (Panel C)
dfKg <- df %>% rename(priceSeries = kgPrice)
# Take the means and standard deviations of the series
dfRiskReward <- processRiskRewardDf(dfKg)
# Drop a few outliers
#dfRiskReward %>% subset(group == "Other") %>% .$Item # "Lupins" "Quinoa" "True hemp, raw or retted" "Vetches"  "Yautia"  
dfRiskReward <- dfRiskReward %>% subset(group != "Other")
# Save for use in Figure 2
dfRRkg <- dfRiskReward; dfRRkg$facetKey <- "C) Producer price (constant 2014-16 USD/kg)" 
# Get slope and y-intercept of the fluctuation scaling and add to panel title
outVec <- getSlopeEtc(dfRiskReward)
m <- outVec[1]; b <- outVec[2]; adjR2 <- outVec[3]; N <- outVec[4];
dfKgPrice <- dfRiskReward
thisFacet <- paste0("N = ", N, " agricultural commodities", "\nAdj. R-squared = ", adjR2, ", Slope = ", m, ", Y intercept = ", b)
dfKgPrice$facetKey <- paste0("C) Producer price (constant 2014-16 USD/kg)\n", thisFacet)
#----------------------------------------------------------------------------
# Subset the data for the yield plot (Panel D)
dfYd <- df %>% rename(priceSeries = Yield)
# Take the means and standard deviations of the series
dfRiskReward <- processRiskRewardDf(dfYd)
# Drop a few outliers
dfRiskReward <- dfRiskReward %>% subset(group != "Other")
# Save for use in Figure 2
dfRRyd <- dfRiskReward; dfRRyd$facetKey <- "D) Yield (kg/ha)"
# Get slope and y-intercept of the fluctuation scaling and add to panel title
outVec <- getSlopeEtc(dfRiskReward)
m <- outVec[1]; b <- outVec[2]; adjR2 <- outVec[3]; N <- outVec[4];
dfYd <- dfRiskReward
thisFacet <- paste0("N = ", N, " agricultural commodities", "\nAdj. R-squared = ", adjR2, ", Slope = ", m, ", Y intercept = ", b)
dfYd$facetKey <- paste0("D) Yield (kg/ha)\n", thisFacet)
#----------------------------------------------------------------------------
# Subset the data for the producer price/hectare plot (Panel B)
dfHa <- df %>% rename(priceSeries = haPrice)
# Take the means and standard deviations of the series
dfRiskReward <- processRiskRewardDf(dfHa)
# Drop a few outliers
dfRiskReward <- dfRiskReward %>% subset(group != "Other")
# Save for use in Figure 2
dfRRha <- dfRiskReward; dfRRha$facetKey <- "B) Producer price (constant 2014-16 USD/ha)"
# Get slope and y-intercept of the fluctuation scaling and add to panel title
outVec <- getSlopeEtc(dfRiskReward)
m <- outVec[1]; b <- outVec[2]; adjR2 <- outVec[3]; N <- outVec[4];
dfHaPrice <- dfRiskReward
thisFacet <- paste0("N = ", N, " agricultural commodities", "\nAdj. R-squared = ", adjR2, ", Slope = ", m, ", Y intercept = ", b)
dfHaPrice$facetKey <- paste0("B) Producer price (constant 2014-16 USD/hectare)\n", thisFacet)
#----------------------------------------------------------------------------
# Collect the processed datasets into a single plot data frame 
dfPlot2 <- do.call(rbind, list(dfKgPrice, dfHaPrice, dfYd)) %>% as.data.frame()
dfPlot2$group <- as.factor(dfPlot2$group)
#----------------------------------------------------------------------------
# Plot panels B-D of Figure 1
names(colorVec) <- levels(dfPlot2$group)
names(shapeVec) <- levels(dfPlot2$group)
fillScale <- scale_fill_manual(name = "group", values = colorVec)
shapeScale <- scale_shape_manual(name = "group", values = shapeVec)
gg <- ggplot(dfPlot2, aes(x = lmPrice, y = lsPrice,
                          group = group,
                          fill = group,
                          shape = group))
gg <- gg + geom_smooth(aes(group = NULL, fill = NULL, shape = NULL), method = lm, se = F)
gg <- gg + geom_point(alpha = 0.6, size = pointSize, color = "black", stroke = 0.5)
gg <- gg + shapeScale
gg <- gg + fillScale
gg <- gg + labs(x = "Nat. log of mean", y = "Nat. log of stand. dev.")
gg <- gg + facet_wrap(~facetKey, nrow = 2, scales = "free")
gg <- gg + theme_bw()
gg <- gg + theme(axis.title = element_text(size = axisTitleSize),
                 axis.text = element_text(size = axisTextSize),
                 legend.position = "bottom", #"right",
                 legend.spacing.x = unit(0.25, 'cm'),
                 legend.text = element_text(size = legendTextSize),
                 legend.title = element_blank(),
                 strip.background = element_rect(fill = "white"),
                 strip.text = element_text(hjust = 0, size = axisTitleSize))
gg <- gg + guides(
  fill = guide_legend(
    ncol = 2,
    byrow = TRUE,
    override.aes = list(linetype = 0, shape = shapeVec)
  ),
  shape = guide_legend(
    ncol = 2
  ),
  color = guide_legend(
    ncol = 2,
    override.aes = list(linetype = 0)
  )
)
gg <- gg + theme(
  legend.position = c(0.8, 0.2), # Adjust coordinates to land in the empty facet
  legend.justification = c(0.6, 0.5)
)
ggFig1BCD <- gg
#----------------------------------------------------------------------------
# Create Figure 1 by putting panels A and B-D together into a single plot
# using the patchwork library
ggFig1 <- (ggFig1A / ggFig1BCD) + 
  plot_annotation(
    title = "Empirical evidence of fluctuation scaling in agricultural commodity prices",
    subtitle = "FAO price and yield series 1961-2024, dollar values in constant 2014-16 USD",
    theme = theme(
      plot.title = element_text(size = titleSize),
      plot.subtitle = element_text(size = subtitleSize)
    )
  )
#===========================================================================
# Show Figure 1
ggFig1
# Use e.g. ggsave("fig1.png", width, height, dpi = 350) to save figure as .png, .tif, etc.
#===========================================================================
# Figure 2, panels A-D - Mean–dispersion relationship in logged price and yield series
#===========================================================================
# Collect the previously saved data frames into single data frame for plotting
dfPlot3 <- do.call(rbind, list(dfRRkg, dfRRha, dfRRyd, dfRRexp)) %>% as.data.frame()
#----------------------------------------------------------------------------
# Plot Fig 2 panels A-D
ggFig2 <- ggplot(dfPlot3, aes(x = mlPrice, y = log(slPrice))) +
  geom_smooth(aes(group = NULL, fill = NULL, shape = NULL), method = lm, se = F) +
  geom_point(size = 0.5) + facet_wrap(~facetKey, scales = "free_x", ncol = 2) +
  labs(x = "Mean of logged series", y = "Standard deviation of logged series (log scale)") +
  theme_bw() +
  theme(axis.title = element_text(size = axisTitleSize),
        axis.text = element_text(size = 6),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(hjust = 0, size = 6))
#===========================================================================
# Show Figure 2
ggFig2
# Use e.g. ggsave("fig1.png", width, height, dpu = 350) to save figure as .png, .tif, etc.
#===========================================================================
# End
#===========================================================================
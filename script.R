#--------------------------------------------------------

suppressPackageStartupMessages(library(smoof))

fn <- function(X){
  if(!is.matrix(X)) X <- matrix(X, nrow = 1) # <- if a single vector is passed as X
  
  Y <- apply(X, MARGIN = 1,
             FUN = smoof::makeRosenbrockFunction(dimensions = dim))
  
  return(Y)
  
}

## Config 1
recpars1 <- list(name = "recombination_exp", cr = 0.6)
mutpars1 <- list(name = "mutation_best", f = 2)

## Config 2
recpars2 <- list(name = "recombination_geo", alpha = 0.6)
mutpars2 <- list(name = "mutation_rand", f = 1.2)

set.seed(2)

dimensions <- sample(2:150, 64, replace=FALSE)
results1 <- matrix(0, 30, 64)

i <- 1
for(d in 1:64) {
  for(t in 1:30){
    dim <- dimensions[d]# d
    selpars <- list(name = "selection_standard")
    stopcrit <- list(names = "stop_maxeval", maxevals = 5000 * dim, maxiter = 100 * dim)
    probpars <- list(name = "fn", xmin = rep(-5, dim), xmax = rep(10, dim))
    popsize = 5 * dim
    
    suppressPackageStartupMessages(library(ExpDE))
    
    # Run algorithm on problem:
    out <- ExpDE(mutpars = mutpars1,
                 recpars = recpars1,
                 popsize = popsize,
                 selpars = selpars,
                 stopcrit = stopcrit,
                 probpars = probpars,
                 showpars = list(show.iters = "dots", showevery = 20))
    
    # Extract observation:
    out$Fbest
    results1[t, d] = out$Fbest
  }
  i = i + 1
}

results2 <- matrix(0, 30, 64)
i <- 1
print(i)
for(d in 1:64) {
  for(t in 1:30){
    dim <- dimensions[d]# d
    
    selpars <- list(name = "selection_standard")
    stopcrit <- list(names = "stop_maxeval", maxevals = 5000 * dim, maxiter = 100 * dim)
    probpars <- list(name = "fn", xmin = rep(-5, dim), xmax = rep(10, dim))
    popsize = 5 * dim
    
    suppressPackageStartupMessages(library(ExpDE))
    
    # Run algorithm on problem:
    out <- ExpDE(mutpars = mutpars2,
                 recpars = recpars2,
                 popsize = popsize,
                 selpars = selpars,
                 stopcrit = stopcrit,
                 probpars = probpars,
                 showpars = list(show.iters = "dots", showevery = 20))
    
    # Extract observation:
    out$Fbest
    
    results2[t, i] = out$Fbest
  }
  i = i + 1
}


#----------------------------------------------------------------------------
# clean workspace
rm(list=ls())

# install required packages if needed
packages_needed <- c("stringr","ggplot2", "multcomp")
for (package_name in packages_needed) {      
  if (!(package_name %in% rownames(installed.packages()))){
    install.packages(package_name)
  }
}
#===================
####
# Load data
####

data <- read.table("data.csv",
                   header = TRUE, sep = ',')

# View a part of the data filtered by group
data[which(data$grou =='1'),][1:12,]
data[which(data$group=='2'),][1:12,]

summary(data)

# Aggregate data (algorithm means by instance group)
aggdata <- with(data,
                aggregate(x   = result,
                          by  = list(group, t),
                          FUN = mean))

# Rename columns
names(aggdata) <- c("Group", 
                    "Instance",
                    "Y")

for (i in 1:2){
  aggdata[, i] <- as.factor(aggdata[, i])
}

# Make factor level names more informative
levels(aggdata$Group) <- c("Config",
                           unlist(lapply("Config",
                                         paste0,
                                         1:2)))

summary(aggdata)

####
## Exploratory data analysis: plot observations by Algorithm and Instance_Group
####

library(ggplot2)

# png(filename = "../figs/algo_lineplot.png",
#     width = 1000, height = 400, 
#     bg = "transparent")
p <- ggplot(aggdata, aes(x = Instance, 
                         y = Y, 
                         group = Group, 
                         colour = Group))
p + geom_line(linetype=2) + geom_point(size=5)
# dev.off()

####
## Statistical modeling
####

# First model
# Response is result of the overall mean plus algorithm effect (treatment)
# and instance (block) effect

model <- aov(Y~Group+Instance,
             data = aggdata)

summary(model)
summary.lm(model)$r.squared

# Graphical test of assumptions
# png(filename = "../figs/algo_res1.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mfrow = c(2, 2))
plot(model, pch = 20, las = 1)
# dev.off()

# Try log-transformed data

model2 <- aov(log(Y)~Group+Instance,
              data = aggdata)
summary(model2)
summary.lm(model2)$r.squared

# Graphical test of assumptions
# png(filename = "../figs/algo_res2.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mfrow = c(2, 2))
plot(model2, pch = 20, las = 1)
# dev.off()

library(car)
# png(filename = "../figs/algo_qq.png",
#     width = 600, height = 600, 
#     bg = "transparent")

# model diagnostic
par(mfrow = c(1, 1))
qqPlot(model2$residuals, pch = 20, las = 1)
# dev.off()

####
## Blocking efficiency
####

mydf        <- as.data.frame(summary(model2)[[1]])
MSblocks    <- mydf["Instance","Mean Sq"]
MSe         <- mydf["Residuals","Mean Sq"]
a           <- length(unique(aggdata$Group))
b           <- length(unique(aggdata$Instance))
((b - 1) * MSblocks + b * (a - 1) * MSe) / ((a * b - 1) * MSe)

####
## Post-hoc multiple comparisons
####

# using model 2 - log-transformed data
library(multcomp)
duntest     <- glht(model2,
                    linfct = mcp(Group = "Dunnett"))

summary(duntest)


duntestCI   <- confint(duntest)
# png(filename = "../figs/algo_mcp.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mar = c(5, 8, 4, 2), las = 1)
plot(duntestCI,
     xlab = "Mean difference (log scale)")
# dev.off()

# using model 1
duntest1     <- glht(model,
                     linfct = mcp(Group = "Dunnett"))

summary(duntest1)

duntestCI1   <- confint(duntest1)

# png(filename = "../figs/algo_mcp.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mar = c(5, 8, 4, 2), las = 1)
plot(duntestCI1,
     xlab = "Mean difference")
# dev.off()

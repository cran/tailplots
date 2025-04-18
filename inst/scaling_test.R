sample_sizes <- c(10, 100, 500)


distributions <- list(
  "Gamma (shape=2, scale=1)" = function(n) rgamma(n, shape = 2, scale = 1),
  "Gamma (shape=0.5, scale=2)" = function(n) rgamma(n, shape = 0.5, scale = 2),
  "Log-normal (meanlog=0, sdlog=1)" = function(n) rlnorm(n, meanlog = 0, sdlog = 1),
  "Weibull (shape=1.5, scale=1)" = function(n) rweibull(n, shape = 1.5, scale = 1),
  "Weibull (shape=0.7, scale=1)" = function(n) rweibull(n, shape = 0.7, scale = 1),
  "Exponential (rate=1)" = function(n) rexp(n, rate = 1),
  "Exponential (rate=0.5)" = function(n) rexp(n, rate = 0.5),
  "Pareto (shape=2, scale=1)" = function(n) 1 / (runif(n)^0.5),
  "Pareto (shape=3, scale=1)" = function(n) 1 / (runif(n)^(1/3)),
  "Uniform (0,5)" = function(n) runif(n, min = 0, max = 5),
  "Uniform (0,10)" = function(n) runif(n, min = 0, max = 10),
  "Beta (shape1=2, shape2=5)" = function(n) rbeta(n, shape1 = 2, shape2 = 5),
  "Beta (shape1=5, shape2=2)" = function(n) rbeta(n, shape1 = 5, shape2 = 2)
)

##
# Gamma
##


# Define the path to the parent directory
parent_directory <- normalizePath("..")  # Moves up one level from the current directory
pdf_path <- file.path(parent_directory, "all_gamma_tailplots.pdf")

# Open a PDF file in the parent directory
pdf(pdf_path)

# Loop over sample sizes and distributions
for (n in sample_sizes) {
  for (dist_name in names(distributions)) {
    x <- distributions[[dist_name]](n)
    
    # Create a new plot
    plot_title <- paste("Distribution:", dist_name, "- Sample Size:", n)
    #"o" functioniert, "b" nicht!
    gamma_tailplot(x, method = "unbiased", xscale = "b")
    #gamma_tailplot(x, method = "unbiased", xscale = "b")
    title(main = plot_title)
  }
}

dev.off()




##
# Pareto
##


# Define the path to the parent directory
parent_directory <- normalizePath("..")  # Moves up one level from the current directory
pdf_path <- file.path(parent_directory, "all_pareto_tailplots.pdf")

# Open a PDF file in the parent directory
pdf(pdf_path)

# Loop over sample sizes and distributions
for (n in sample_sizes) {
  for (dist_name in names(distributions)) {
    x <- distributions[[dist_name]](n)
    
    # Create a new plot
    plot_title <- paste("Distribution:", dist_name, "- Sample Size:", n)
    #"o" functioniert, "b" nicht!
    pareto_tailplot(x, method = "unbiased", xscale = "b")
    #pareto_tailplot(x, method = "unbiased", xscale = "b")
    title(main = plot_title)
  }
}

dev.off()

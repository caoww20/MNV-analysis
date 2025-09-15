

# Extract AAAS color palette
color1 = pal_aaas()(10)
show_col(color1)
# Nature-style color palette
colorPalette<-c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF")
# Heatmap color palette
my_colors <- colorRampPalette(c("#9AC9E0","#2272B4","#1A325F","#040419"))(100)

my_colors <- colorRampPalette(colors = c("#F06143","#E33641","#AE1758","#471D49","#040419"))(100)
# PCA plot color palette
my_colors <- c("#7E277B", "#FCA41C","#0D8040","#3353A3","#EF251F")
# Custom color palettes
## Gray
colorPalette<-c('#F4F4F4','#DBDBDB','#A0A0A0','#606060','#424242')
## Blue tones
colorPalette<-c( "black","#1A325F","#2A5888","#4396B1","#89CEED","#C4F5FC")
## Bar plot (ignore for now)
# colorPalette<-c( "#1A325F","#2E6AA4","#55ADD2","#82C0C5","#BEE1C6")
## Red and blue
colorPalette<-c("#B0282E","#2A5888")
## Other schemes
colorPalette<-c("#397F43","#4daf4a")

colorPalette<-c('#B0282E',"#D03650","#e41a1c","#F0B598")

# Set random seed for reproducibility (optional)

set.seed(123)

# Generate a random matrix
df <- matrix(runif(2*5000), nrow = 5000, ncol = 2)
df <- as.data.frame(df)
colnames(df)<-c('x1','x2')
df$x1x2<-df$x1*df$x2
df$y<-sample(0:1, 5000, replace = TRUE)

# Measure execution time
{
  time_taken <- system.time({
    # Fit a logistic regression model
    model <- glm(y ~ x1 + x2 + x1x2, data = df, family = binomial)
  })[3]
  # Print model summary
  summary(model)
  print(time_taken)
}

# Measure execution time
{
  
  time_taken <- system.time({
    # Create a sample vector x (replace with your own)
    model_list=list()
    x <- 1:5000
    for (i in 1:10) {
      
      # Randomly sample 50 elements from x
      random_sample <- sample(x, 50)
      
      # Fit a logistic regression model
      model <- glm(y ~ x1 + x2 + x1x2, data = df[random_sample,], family = binomial)
      
      # Save model summary
      model_list[[i]]=summary(model)
    }
  })[3]
  # Print elapsed time
  print(time_taken)
}


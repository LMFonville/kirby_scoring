#### KIRBY SCORING TOOL ####


#### Prep Reference Dataframe ####

# Amount 1 represents first option on questionnaire and amount 2 is the second option
# Same applies to delay1 and delay2
itemdf <- data.frame(
  item = seq(1,27,1),
  amount1=c(55,55,19,85,25,50,15,60,78,40,30,67,35,50,69,49,85,24,80,28,34,30,41,60,54,22,55),
  amount2=c(54,75,25,31,14,47,35,25,80,55,11,75,34,27,85,60,80,35,33,30,50,25,75,54,80,25,20),
  delay1=c(117,0,0,7,19,160,0,14,0,0,7,0,186,21,0,0,157,0,14,0,0,80,0,111,0,0,7),
  delay2=c(0,61,53,0,0,0,13,0,162,62,0,119,0,0,91,89,0,29,0,179,30,0,20,0,30,136,0),
  # calculate k and size of magnitude below
  k = rep(0, 27), size = NA)

# A is the larger, delayed amount
# V is the smaller immediate reward
# D is the delay associated with the larger amount A
# Size is the magnitude of the larger amount A
# Small: A<=£35 (options are 25,30,35)
# Medium: A>£35 and A<=£60 (options are 50,55,60)
# Large: A>£60 and A<=£85 (options are 75,80,85)

# Use for loop to stick with base R
for (i in 1:nrow(itemdf)){
  # Calculate degree of discounting (k) at 'indifference' for each item
  # These values match what I found online for fixed values of k
  itemdf$k[i] <- ((max(itemdf$amount1[i], itemdf$amount2[i]) / min(itemdf$amount1[i], itemdf$amount2[i])) - 1) / 
    max(itemdf$delay1[i], itemdf$delay2[i])
  # Specify size of delayed reward as small, medium, large 
  # These match what I found online as well
  m = max(itemdf$amount1[i], itemdf$amount2[i])
  if (m<=35){itemdf$size[i]=1} else if(m>35 & m <=60){itemdf$size[i]=2} else if(m>60 & m <=85){itemdf$size[i]=3}
  # calculate reward ratio as well
  itemdf$rr[i] <- min(itemdf$amount1[i], itemdf$amount2[i]) / max(itemdf$amount1[i], itemdf$amount2[i])
}
itemdf$size <- factor(itemdf$size, levels = c(1, 2, 3), labels = c('S', 'M', 'L'))
# a roughly similar value of k is provided for each magnitude of delayed reward
# order the data based on k and the magnitude of the reward in a fixed order of small, medium, large
itemdf$order <- 0
# assign rank to increasing values of k within size and add margin to keep the order consistent
itemdf$order[itemdf$size=='S'] <- rank(itemdf$k[itemdf$size=='S'])+0.1
itemdf$order[itemdf$size=='M'] <- rank(itemdf$k[itemdf$size=='M'])+0.2
itemdf$order[itemdf$size=='L'] <- rank(itemdf$k[itemdf$size=='L'])+0.3

#### Specify functions for calculations ####
# function to calculate values of k as well as consistency
estimate.k <- function(df){
  # recode participant response into choice
  # if response is 1 (amount1) and duration (delay1) is zero, or if response is 2 (amount2) and duration (delay2) is zero, the immediate reward was selected
  for (i in 1:nrow(df)){
    if (df$resp[i]==1 & df$delay1[i]==0){
      df$choice[i] <- 0 # smaller immediate reward response
    } else if (df$resp[i]==2 & df$delay2[i]==0){
      df$choice[i] <- 0 # smaller immediate reward response
    } else {
      df$choice[i] <- 1 # larger delayed reward response
    }
  }
  # at smallest k, count how many times delayed reward is selected
  # at each increasing k-order, count how many previous items the short reward was chosen 
  # and for how many in the current and subsequent the delayed reward was chosen
  #k_df <- df[,c('item', 'k', 'choice')]
  df_k <- df
  df_k$k_count = 0
  # first index counts number of times delayed response was chosen
  df_k$k_count[1] = sum(df_k$choice)
  for (i in 2:nrow(df_k)){
    # each subsequent item considers:
    # the number of times the immediate reward was chosen prior to current value of k
    # and number of times the delayed reward was chosen for subsequent values of k
    prior=sum(df_k$choice[1:(i-1)]==0)
    remaining=sum(df_k$choice[i:nrow(df_k)])
    df_k$k_count[i] = prior + remaining
  }
  # append with row to consider if all choices were for immediate reward
  df_k[nrow(df_k)+1,] <- NA
  df_k$k[nrow(df_k)] <- max(df_k$k, na.rm = TRUE)
  df_k$k_count[nrow(df_k)] <- sum(df_k$choice==0, na.rm=TRUE)
  # find highest consistency of choice (when most previous choices are immediate reward and subsequent choices are delayed reward)
  consistency <- max(df_k$k_count) 
  # find indices where consistency was highest
  c_idx <- which(df_k$k_count==consistency)
  # if index contains 1 (all choices for delayed reward), select lowest value
  # if index contains maximum row number (all choices for immediate reward), select largest value
  geo_k <- NULL
  for (j in c_idx){
    geo_k <- c(geo_k, ifelse(j==1, df_k$k[j], sqrt(df_k$k[j]*df_k$k[j-1])))
  }
  if (length(geo_k)>1){
    geo_k <- prod(geo_k)^(1/length(c_idx))
  }
  # calculate proportion of choices
  all1 <- sum(df$choice==0)/nrow(df)
  all2 <- sum(df$choice)/nrow(df)
  
  # produce logistic regression output as well (not fully tested!)
  # calculate the probability of choosing the delayed reward (response 1) using reward ratio and delay time
  # reward ratio is immediate amount relative to delayed amount (calculated in the initial loop)
  # regressor coefficients are 1/(1-R) and delay time
  df$r1 <- 1-(1/df$rr)
  #df$r1 <- 1-df$rr
  df$r2 <- apply(df[4:5], 1, function(x) max(x))
  k.glm <- glm(choice ~ 0 + r1 + r2, data = df, family = binomial)
  k.model <- as.numeric(k.glm$coefficients[2] / k.glm$coefficients[1])
  # Include McFadden's R-squared value
  rsquared <- with(summary(k.glm), 1 - deviance/null.deviance)
  # collate relevant outputs
  k.output <- list(k.response=geo_k, consistency=consistency, 
                   consistency_pct = consistency / nrow(df),
                   choiceNow=all1, choiceDelay=all2,
                   glm=k.glm, k.model=k.model, model.r2=rsquared, data_k=df_k, data_glm=df)
  return(k.output)
}

df.to.k <- function(itemdf, sub.df){
  # only run with complete data, otherwise skip
  subj=unique(sub.df$participant)
  if (nrow(sub.df)<27){
    return(paste0("Incomplete data for sub-NK", subj))
  }
  # combine itemdf with subject scores
  sub.itemdf <- itemdf
  sub.itemdf$resp <- sub.df$score
  # reorder data according to increasing values of k and fixed sequence for magnitude of delayed reward
  sub.itemdf <- sub.itemdf[order(sub.itemdf$order),]
  # calculate relevant metrics
  sub.kAll <- estimate.k(sub.itemdf)
  sub.kSmall <- estimate.k(sub.itemdf[sub.itemdf$size=='S',])
  sub.kMedium <- estimate.k(sub.itemdf[sub.itemdf$size=='M',])
  sub.kLarge <- estimate.k(sub.itemdf[sub.itemdf$size=='L',])
  kCombined <- prod(sub.kSmall$k.response, sub.kMedium$k.response, sub.kLarge$k.response)^(1/3)
  # collate the values of k (overall, small, medium, large) and their consistency for now
  # there's other variables that are calculated but unlikely to be important
  sub.k.df <- data.frame(id=subj, kOverall=sub.kAll$k.response, kCombined=kCombined, kSmall=sub.kSmall$k.response,
                      kMedium=sub.kMedium$k.response,kLarge=sub.kLarge$k.response,
                      consistencyOverall=sub.kAll$consistency, consistencySmall=sub.kSmall$consistency, 
                      consistencyMedium=sub.kMedium$consistency,consistencyLarge=sub.kLarge$consistency,
                      proportionOverall=mean(sub.kAll$data_k$choice, na.rm = TRUE), proportionSmall=mean(sub.kSmall$data_k$choice, na.rm = TRUE),
                      proportionMedium=mean(sub.kMedium$data_k$choice, na.rm = TRUE), proportionLarge=mean(sub.kLarge$data_k$choice, na.rm = TRUE))
  return(sub.k.df)
}

#### Single file usage ####
# adjust below to work with input data
# make sure input df columns match itemdf
sub.file <- file.choose()
sub.df <- read.csv(sub.file)
# warning messages to do with glm.fit can be ignored as we don't use those outputs
sub.k <- df.to.k(itemdf, sub.df)

#### Loop through files in a folder ####

# this can easily be applied in a loop as well
input.dir <- choose.dir()
# change pattern to match data
KDD_files <- list.files(path=input.dir, pattern='KDD', recursive = FALSE, full.names = TRUE) # don't use full names to get easy access to subject ID in filename
# maybe add a check to ensure only the right files are found in KDD_files?

# loop through files and get relevant information
kirbydf <- NULL
for (subfile in KDD_files){
  sub.df <- read.csv(subfile)
  # warning messages to do with glm.fit can be ignored as we don't use those outputs
  sub.k <- df.to.k(itemdf, sub.df)
  # simply bind rows to produce single dataframe
  kirbydf <- rbind(kirbydf, sub.k)
}


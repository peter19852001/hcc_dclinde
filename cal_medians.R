
# to calculate the medians of the columns of a csv file

# assume infile and outfile have been defined
cal.medians <- function(data) {
  # data is a dataframe
  # out is output file name
  n <- ncol(data);
  r <- rep(0,n);
  for(i in 1:n) {
    r[i] <- median(data[,i], na.rm=TRUE);
  }
  r
}

tmp.data <- read.csv(file=infile, header=TRUE);
tmp.medians <- cal.medians(tmp.data);
cat(tmp.medians, file=outfile, sep=",");
cat("\n");

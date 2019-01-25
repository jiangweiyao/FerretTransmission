#############################
# Initial set-up of the R session
rm(list=ls())
options(stringsAsFactors = F)

##################################
# Specify input file and read it
ferret.dir="D:/SPounds/BioResearch/InfectiousDisease/Rosch/FerretModel/"
ferret.file=paste(ferret.dir,"ferret read matrix.csv",sep="")
ferret.data0=read.csv(ferret.file,header=T,as.is=T)

####################################
# Specify name of result files
res.dir="D:/SPounds/BioResearch/InfectiousDisease/Rosch/Transmission-Manuscript-2018-07-13/Resubmit-2018-10/"
stat.result.csv.file=paste0(res.dir,"/ferret-read-count-stats.csv")
plot.file=paste0(res.dir,"/ferret-read-count-plots.pdf")

##################################
# Copy input file to result directory
file.copy(ferret.file,res.dir)

#################################
# standardize column names
colnames(ferret.data0)=gsub(".","",colnames(ferret.data0),fixed=T)
colnames(ferret.data0)=gsub("recipeint","recipient",colnames(ferret.data0),fixed=T)

#################################
# Generate plots to justify choice of read count = 10 for defining infection

# Identify columns with read count data
input.index=grep("input",colnames(ferret.data0))
donor.index=grep("donor",colnames(ferret.data0))
rcpt.index=grep("recipient",colnames(ferret.data0))

# Find rows with zero input
input.sum=rowSums(ferret.data0[,input.index])   # compute input total for each OTU (matrix row)
zero.input=which(is.element(input.sum,c(0,NA))) # identify rows with 0 or NA as the input total

# open plot file for capturing plots
pdf(plot.file,height=10,width=8)
###########################################################
# Generate plots of donor counts for rows with zero input
par(mfrow=c(3,3))
for (i in donor.index) # loop over donors
{
  donor.number=gsub("donor","",colnames(ferret.data0)[i]) # identify column with data for this donor
  donor.count=table(ferret.data0[zero.input,i])           # number of OTUs with specific number of reads
  x=as.numeric(names(donor.count))                        # read count for OTUs
  plot(c(0,10),c(0,length(zero.input)),type="n",          # generate plot space
       xlab="Donor Read Count",
       ylab="Number of OTUs",
       main=paste0("Donor ",donor.number,
                   " Counts for OTUs with 0 Input"))      # add histogram
  lines(x,donor.count,type="h",lwd=3)
} # end loop over donors


###############################################################
# Generate plots of number of OTUs with a given read count for all donors
for (i in donor.index) # loop over donors
{
  donor.number=gsub("donor","",colnames(ferret.data0)[i]) # identify column with data for this donor
  donor.count=table(ferret.data0[,i])                     # tabulate number of OTUs with each read count
  y=log10(as.numeric(names(donor.count))+1)               # log10 transform for plot
  x=as.vector((donor.count))                              # read count values
  plot(c(0,50),c(0,7),xlab="Read Count",                  # generat plot space
       ylab="log10(Number of OTUs+1)",
       type="n",main=paste0("Donor ",donor.number))
  rect(0,0,10.5,7,col="gainsboro",border=NA)              # add background rectangle
  lines(x,y,type="h",lwd=3)                               # add histogram
} # end loop over donors


###############################################################
# Generate plots of number of OTUs with a given read count for all recipients
for (i in rcpt.index) # loop over recipients
{
  rcpt.number=gsub("recipient","",colnames(ferret.data0)[i]) # identify column with data for this recipient
  rcpt.count=table(ferret.data0[,i])                         # tabulate data for this recipient
  y=log10(as.numeric(names(rcpt.count))+1)                   # log10 transform for plot
  x=as.vector(rcpt.count)                                    # read count values
  plot(c(0,100),c(0,7),xlab="Read Count",                    # generate plot space
       ylab="log10(Number of OTUs+1)",
       type="n",main=paste0("Recipient ",rcpt.number))
  rect(0,0,10.5,7,col="gainsboro",border=NA)                 # add gray background rectangle
  lines(x,y,type="h",lwd=3)                                  # add histogram
} # end loop over recipients

dev.off() # close plot file


###############################################################
# Compute result table
number.donors.infected=rowSums(ferret.data0[,donor.index]>=10) # number of donors infected with each OTU

# mean probability of transmission
mean.prob.transmit=rep(0,length(number.donors.infected)) # initialize result   
for (i in donor.index) # loop over donors
{
  donor.number=gsub("donor","",colnames(ferret.data0)[i])                 # identify column with data for this donor
  rcpt.indx=grep(paste0("recipient",donor.number),colnames(ferret.data0)) # identify recipients of this donor
  pr.trans=rowMeans(ferret.data0[,rcpt.indx]>10)                          # compute mean transmission probability for this donor's recipients
  pr.trans=pr.trans*(ferret.data0[,i]>10)                                 # donor must have infection for transmission to occur
  mean.prob.transmit=mean.prob.transmit+pr.trans                          # add that to the total
} # end loop over donors

mean.prob.transmit=mean.prob.transmit/number.donors.infected # divide total by number of infected donors
mean.prob.transmit[number.donors.infected==0]=0              # set transmission to zero for OTUs with no infected donors

# annotate results in a data.frame
final.result=cbind.data.frame(ferret.data0[,c("locusID","annotation")],
                              number.donors.infected=number.donors.infected,
                              mean.prob.recipient.infected=mean.prob.transmit)

# order results by probability of recipient infection and number of donors infected
result.ord=order(final.result$mean.prob.recipient.infected,
                 -final.result$number.donors.infected)

final.result=final.result[result.ord,]

# write result table file
write.csv(final.result,stat.result.csv.file,
          row.names=F,quote=T)

#!/usr/bin/Rscript
library(reshape2)
library("dplyr")
infile="/rgs01/project_space/roschgrp/TnSeq/partners/CMPB/WorkSpace/TNseq_Hannah5_ferret/s07_tnseq_analysis_BHN97/tnseq.total.hits.tsv"
sam="/rgs01/project_space/roschgrp/TnSeq/partners/CMPB/WorkSpace/TNseq_Hannah5_ferret/pheno.tsv"
gff="/rgs01/project_space/roschgrp/SNP_JR/partners/CMPB/WorkSpace/SNP_BHN97/BHN97.gff3"


### functions ###
ann_table <- function(df, ann, fld1="pos", fld2="pos"){
	getValue <- function(x, data, rng_start="Start", rng_end="End", return_field="Comment"){
		### check if value overlapped and return the first match
		tmp <- data.frame(data %>%
			filter(Start <= x, x <= End) %>%
			filter(row_number() == 1) %>%
			mutate(pos=x)
		)
		if(nrow(tmp) > 0){
			return(tmp)
		}

	}
	
	#df[[fld]] <- as.numeric(gsub("^.*\\.","", rownames(df)))
	df[[fld1]] <- as.numeric(df[[fld1]])
	res <- lapply(df[[fld1]], getValue, data=ann) %>% bind_rows()
	### merge table
	total <- unique(merge(df, res, by.x=fld1, by.y=fld2, all.x=TRUE))
	return(total)
	
}
##################


data=read.table(infile,sep="\t",header=T)

### load sample data
samples <- read.table(sam, sep="\t", header=T, stringsAsFactors=FALSE) ### 63, 8
samples_ex_room_input <- samples[-grep("room|input", samples$ferret),] ### 52, 8


### merge tables
data_cast <- dcast(data, Genome+Position~Sample, fun.aggregate = mean, value.var='Counts')
data_cast[is.na(data_cast)] <- 0
write.table(data_cast, file='tnseq.total.hits.cast.table.tsv',sep='\t', quote=F,row.name=F)



### to keep those samples in the sample sheet
data_castf <- data_cast[,which(colnames(data_cast) %in% c('Genome', 'Position', as.character(samples_ex_room_input$sample.name)))]
write.table(data_cast, file='tnseq.total.hits.cast.table.only.compared.tsv',sep='\t', quote=F,row.name=F)
rownames(data_castf) <- paste(data_castf$Genome, data_castf$Position)

dataf <- data[which(data$Sample %in% samples_ex_room_input$sample.name),]

### add pheno to dataf
dataf_pheno <- merge(dataf, samples, by.x="Sample", by.y="sample.name", all.x=TRUE) ### 24761, 11


### pool
ferret_pool <- data.frame(dataf_pheno %>% group_by(group, ferret, donor.recipient, Genome, Position) %>%
	summarise(pooledCounts=sum(Counts)))

ferret_pool.group.summary <- data.frame(ferret_pool %>% 
	group_by(group, Genome, Position) %>%
	summarise(donor_sample_count=sum(donor.recipient=="D"),
	          donor_total_read_count=sum(pooledCounts[donor.recipient=="D"]),
			  recipient_sample_count=sum(donor.recipient=="R"),
			  recipient_total_read_count=sum(pooledCounts[donor.recipient=="R"])
			  )
	) %>%
	mutate(sharing_status=ifelse(donor_sample_count > 0 & recipient_sample_count ==0, "D",
	                      ifelse(donor_sample_count == 0 & recipient_sample_count >0, "R",
	                      "share")
	))
write.table(ferret_pool.group.summary, file='ferret_pool.group.summary.tsv',sep='\t', quote=F,row.name=F)


group_pool <- data.frame(dataf_pheno %>% group_by(group, donor.recipient, Genome, Position) %>%
	summarise(pooledCounts=sum(Counts)))

group_pool.summary <- data.frame(group_pool %>% 
	group_by(Genome, Position) %>%
	summarise(donor_sample_count=sum(donor.recipient=="D"),
	          donor_total_read_count=sum(pooledCounts[donor.recipient=="D"]),
			  recipient_sample_count=sum(donor.recipient=="R"),
			  recipient_total_read_count=sum(pooledCounts[donor.recipient=="R"])
			  )
	) %>%
	mutate(sharing_status=ifelse(donor_sample_count > 0 & recipient_sample_count ==0, "D",
	                      ifelse(donor_sample_count == 0 & recipient_sample_count >0, "R",
	                      "share")
	))
write.table(group_pool.summary, file='group_pool.summary.tsv',sep='\t', quote=F,row.name=F)


### summary
ferret_pool.group.summary.count <- data.frame(ferret_pool.group.summary %>%
	group_by(group) %>%
	summarise(total_sites=n(),
			  donor_only_sites_num=sum(sharing_status=="D"),
			  recipient_only_sites_num=sum(sharing_status=="R"),
			  shared_sites_num=sum(sharing_status=="share")
			  )
	)
rowSums(ferret_pool.group.summary.count[,3:5])

write.table(ferret_pool.group.summary.count, file='ferret_pool.group.summary.count.tsv',sep='\t', quote=F,row.name=F)

group_pool.summary.count <- data.frame(group_pool.summary %>%
	summarise(total_sites=n(),
			  donor_only_sites_num=sum(sharing_status=="D"),
			  recipient_only_sites_num=sum(sharing_status=="R"),
			  shared_sites_num=sum(sharing_status=="share")
			  )
	)
rowSums(group_pool.summary.count[,2:4])

write.table(group_pool.summary.count, file='group_pool.summary.count.tsv',sep='\t', quote=F,row.name=F)


### annotate tables
anndata=read.table(gff,header=FALSE, sep="\t", quote="");
colnames(anndata) <- c("Chr", "Source", "FakeFeature", "Start", "End", "Score", "Strand", "Phase","Comment");
anndata$SimpleComment <- gsub(".*?(ID=.*?);.*?(product=.*?;).*","\\1;\\2",anndata$Comment,perl=TRUE)

ferret_pool.group.summary.ann <- unique(ann_table(ferret_pool.group.summary, anndata, fld1="Position"))
write.table(ferret_pool.group.summary.ann, file='ferret_pool.group.summary.ann',sep='\t', quote=F,row.name=F)


group_pool.summary.ann <- ann_table(group_pool.summary, anndata, fld1="Position")
write.table(group_pool.summary.ann, file='group_pool.summary.ann',sep='\t', quote=F,row.name=F)


### annotate table
data_cast_ann <- ann_table(data_cast, anndata, fld1="Position")
data_cast_ann_gene <- data.frame(data_cast_ann %>%
	group_by(Comment, Start, End, Strand, Phase) %>%
	summarise_each(funs(sum), matches("^1|^input|^room"))
	)
write.table(data_cast_ann_gene, file='data_cast_ann_gene', sep='\t', quote=F,row.name=F)


ferret_pool$group_collapse <- paste(ferret_pool$group, ferret_pool$donor.recipient, sep=".")
### Note after annotation the dataframe may contain duplicated rows
ferret_pool_ann <- unqiue(ann_table(ferret_pool, anndata, fld1="Position"))
ferret_pool_ann_gene <- data.frame(ferret_pool_ann %>%
	group_by(Comment, Start, End, Strand, group_collapse) %>%
	summarise(pooledCounts=sum(pooledCounts))
	)
ferret_pool_ann_gene_cast <- dcast(ferret_pool_ann_gene, Comment+Start+End+Strand~group_collapse, value.var="pooledCounts")
ferret_pool_ann_gene_cast[is.na(ferret_pool_ann_gene_cast)] <- 0
write.table(ferret_pool_ann_gene_cast, file='ferret_pool_ann_gene_cast', sep='\t', quote=F,row.name=F)








############# archived codes
### analysis 1
### group1
#group='G1'
#group='G2'
group='G3'
samples_g <- samples[samples$group == group,]
g <- samples_g$samples
data_castf_G <- data_castf[,which(colnames(data_castf) %in% g)]
data_castf_G <- data_castf_G[which(rowSums(data_castf_G) != 0),]
unique_grp <- unique(samples_g$group2)
add_d_cols <- ''
add_r_cols <- c()
for ( i in 1:length(unique_grp)){
	grps=unique_grp[i]
	print(grps)
	grp_sample=samples_g$samples[which(samples_g$group2 %in% grps)]
	print(grp_sample)
	coln_for_add=paste(grps, grp_sample[1], grp_sample[2], sep='.')
	data_castf_G[coln_for_add]=data_castf_G[grp_sample[1]]+data_castf_G[grp_sample[2]]
	if (grepl('_R_', coln_for_add)){
		add_r_cols[length(add_r_cols)+1]=coln_for_add
	} else {
		add_d_cols=coln_for_add
	}		
}
head(data_castf_G)
for ( i in 1: length(add_r_cols)){
	#comparison=paste(add_r_cols[i], 'compare', sep='_')
	print(add_r_cols[i])
	print(add_d_cols)
	comparison=gsub("\\..*","", add_r_cols[i])	
	comparison=paste(comparison,'compare', sep='_')
	print(comparison)
	data_castf_G[comparison]=ifelse((data_castf_G[add_r_cols[i]] > 0 & data_castf_G[add_d_cols] > 0) ,"RD",ifelse((data_castf_G[add_r_cols[i]] == 0 & data_castf_G[add_d_cols] > 0),"D", ifelse((data_castf_G[add_r_cols[i]] > 0 & data_castf_G[add_d_cols] == 0),"R","-")))
}
write.table(data_castf_G, file=paste('tnseq.total.hits.cast.table.only.compared', group, 'tsv', sep='.'),sep='\t', quote=F,row.name=T)

### analysis 2
doner_ls <- samples$samples[samples$donor_recipient=='D']
data_castf_G <- data_castf
data_castf_G$Genome <- NULL
data_castf_G$Position <- NULL
data_castf_G <- data_castf_G[which(rowSums(data_castf_G) != 0),]
add_d_cols=paste('Doner', paste(doner_ls, collapse='.'), sep='.')
data_castf_G[add_d_cols] <- rowSums(data_castf_G[,colnames(data_castf_G) %in% doner_ls])
unique_grp <- unique(samples$group2[which(samples$donor_recipient=='R')])
for ( i in 1:length(unique_grp)){
        grps=unique_grp[i]
        print(grps)
        grp_sample=samples$samples[which(samples$group2 %in% grps)]
        print(grp_sample)
        coln_for_add=paste(grps, grp_sample[1], grp_sample[2], sep='.')
        data_castf_G[coln_for_add]=data_castf_G[grp_sample[1]]+data_castf_G[grp_sample[2]]
        if (grepl('_R_', coln_for_add)){
                add_r_cols[length(add_r_cols)+1]=coln_for_add
        } else {
                add_d_cols=coln_for_add
        }
}
for ( i in 1: length(add_r_cols)){
        #comparison=paste(add_r_cols[i], 'compare', sep='_')
        print(add_r_cols[i])
        print(add_d_cols)
        comparison=gsub("\\..*","", add_r_cols[i])
        comparison=paste(comparison,'compare', sep='_')
        print(comparison)
        data_castf_G[comparison]=ifelse((data_castf_G[add_r_cols[i]] > 0 & data_castf_G[add_d_cols] > 0) ,"RD",ifelse((data_castf_G[add_r_cols[i]] == 0 & data_castf_G[add_d_cols] > 0),"D", ifelse((data_castf_G[add_r_cols[i]] > 0 & data_castf_G[add_d_cols] == 0),"R","-")))
}
write.table(data_castf_G, file=paste('tnseq.total.hits.cast.table.only.compared', 'question2', 'tsv', sep='.'),sep='\t', quote=F,row.name=T)


### analysis 3
doner_ls <- samples$samples[samples$donor_recipient=='D']
recipient_ls <- samples$samples[samples$donor_recipient=='R']
data_castf_G <- data_castf
data_castf_G$Genome <- NULL
data_castf_G$Position <- NULL
data_castf_G <- data_castf_G[which(rowSums(data_castf_G) != 0),]
add_d_cols=paste('Doner', paste(doner_ls, collapse='.'), sep='.')
add_r_cols=paste('Recipient', paste(recipient_ls, collapse='.'), sep='.')
data_castf_G[add_d_cols] <- rowSums(data_castf_G[,colnames(data_castf_G) %in% doner_ls])
data_castf_G[add_r_cols] <- rowSums(data_castf_G[,colnames(data_castf_G) %in% recipient_ls])
data_castf_G$comparison=ifelse((data_castf_G[add_r_cols] > 0 & data_castf_G[add_d_cols] > 0) ,"RD",ifelse((data_castf_G[add_r_cols] == 0 & data_castf_G[add_d_cols] > 0),"D", ifelse((data_castf_G[add_r_cols] > 0 & data_castf_G[add_d_cols] == 0),"R","-")))
data_castf_G$position <- as.numeric(gsub(".* ","", rownames(data_castf_G)))
data_castf_G$genome <- gsub(" .*","", rownames(data_castf_G))
data_castf_G$gene <- '.'

# read in the gff3
ann='/rgs01/project_space/roschgrp/SNP_JR/partners/CMPB/WorkSpace/SNP_BHN97/BHN97.gff3'
ann_df=read.table(ann, sep='\t', stringsAsFactors=FALSE, fill=T, quote="")
ann_search_start=1
for (i in 1:dim(data_castf_G)[1]){
	print(data_castf_G$position[i])
	insert_pos=data_castf_G$position[i]
	chr=data_castf_G$genome[i]
	for (j in ann_search_start:dim(ann_df)[1]){
		if ( ann_df$V4[j] <= insert_pos & insert_pos <= ann_df$V5[j] & chr == ann_df$V1[j]){
			data_castf_G$gene[i]=ann_df$V9[j]
			ann_search_start=j
			next	
		}
	}
}
write.table(data_castf_G, file=paste('tnseq.total.hits.cast.table.only.compared', 'question3', 'tsv', sep='.'),sep='\t', quote=F,row.name=T)

library(dplyr)
data_castf_G_gene <- data.frame(data_castf_G %>% group_by(gene) %>%
	summarise_(.dots=list(setNames(paste0('sum(', add_d_cols,')'), 'Doner'), setNames(paste0('sum(', add_r_cols,')'), 'Recipient'))))
colnames(data_castf_G_gene) <- c('gene', 'doner', 'recipient')
	#summarise(Doner=sum(Doner.HMR_032.HMR_041.HMR_036.HMR_042.HMR_040.HMR_043), Recipient=sum(Recipient.HMR_029.HMR_055.HMR_030.HMR_056.HMR_031.HMR_057.HMR_033.HMR_058.HMR_034.HMR_059.HMR_035.HMR_060.HMR_037.HMR_061.HMR_038.HMR_062.HMR_039.HMR_063))
write.table(data_castf_G_gene, file=paste('tnseq.total.hits.cast.table.only.compared', 'question3_gene', 'tsv', sep='.'),sep='\t', quote=F,row.name=F)



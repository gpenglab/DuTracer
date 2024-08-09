##define function to handle multiple samples
shannon = function (freq) {
    sum_freq = sum(freq)
    shan = 0
    for (i in freq) {
        index = -(i/sum_freq)*log2(i/sum_freq)
        shan = shan+index
    }
    return(shan)
}


#define a function for describing how many cells share the same intBCs
intBC_across_cell = function(umi_table, color, ylimit=2000, sample_name = "") {
    intBC_cell = umi_table %>% group_by(intBC) %>% summarise(intBCcell = n_distinct(cellBC))
    print(summary(intBC_cell$intBCcell))
    meanBC = mean(intBC_cell$intBCcell)
    intBC_cell = as.data.frame(table(intBC_cell$intBCcell))
    colnames(intBC_cell) = c("intBC_across_cell","freq" )
    nrows = nrow(intBC_cell)
    intBC_cell = intBC_cell[c(1:8, (nrows-1):nrows),] #select 10 important rows
    ggplot(intBC_cell, aes(intBC_across_cell, freq))+geom_bar(fill = color, stat="identity")+theme_classic()+
       theme(axis.text.x = element_text(angle = 90))+ggtitle(sample_name)+xlab("")+
       #geom_text(label = paste0("mean ", round(meanBC, 2)), x=meanBC+2, y=0.7*ylimit)+
       scale_y_continuous(expand=c(0,0),limits = c(0,ylimit))+
       #geom_vline(xintercept=meanBC, lty=2)+
       theme(axis.text = element_text(size=11, color="black"),axis.title=element_text(size=15,color="black"))
}


#define a function for describing how many intBCs integrated in one cell.
intBC_per_cell = function(umi_table, color, ylimit=1000, sample_name = "") {
    intBC_number = umi_table %>% group_by(cellBC) %>% summarise(intBCnumber = n_distinct(intBC)) 
    print(summary(intBC_number$intBCnumber))
    medianBC = median(intBC_number$intBCnumber)
    intBC_number = as.data.frame(table(intBC_number$intBCnumber))
    colnames(intBC_number) = c("intBC_per_cell","freq" )
    selected_nrows = ifelse(nrow(intBC_number)>20, 20, nrow(intBC_number))
    intBC_number = intBC_number[c(1:selected_nrows),] #select 20 important rows
    ggplot(intBC_number, aes(intBC_per_cell, freq))+geom_bar(fill =color, stat="identity")+theme_classic()+
      theme(axis.text.x = element_text(angle = 90))+ggtitle(sample_name)+xlab("")+
     geom_text(label = paste0("median ", round(medianBC, 2)), x=medianBC+2, y=0.7*ylimit)+
      scale_y_continuous(expand=c(0,0),limits = c(0,ylimit))+geom_vline(xintercept=medianBC, lty=2)+
      theme(axis.text = element_text(size=11, color="black"),axis.title=element_text(size=15,color="black"))
}


#define a function for describing how many intBCs integrated in one cell.
intBC_per_cell_violin = function(umi_table, color, ylimit=1000, sample_name = "") {
    intBC_number = umi_table %>% group_by(cellBC) %>% summarise(intBCnumber = n_distinct(intBC)) 
    print(summary(intBC_number$intBCnumber))
    #medianBC = median(intBC_number$intBCnumber)
    #intBC_number = as.data.frame(table(intBC_number$intBCnumber))
    #colnames(intBC_number) = c("intBC_per_cell","freq" )
    intBC_number$sample = sample_name
    #selected_nrows = ifelse(nrow(intBC_number)>column, column, nrow(intBC_number))
    #intBC_number = intBC_number[c(1:selected_nrows),] #select 20 important rows
    p = ggplot(intBC_number, aes(sample, intBCnumber))+geom_violin(fill =color)+theme_classic()+
      theme(axis.text.x = element_text(angle = 90))+ggtitle(sample_name)+xlab("")+ylab("intBCs_per_cell")+
      scale_y_continuous(expand=c(0,0),limits = c(0,ylimit))+
      geom_boxplot(width=0.3)+
      #geom_text( y=0.95*ylimit, label=paste0("medianBC=",medianBC))+
      #geom_hline(yintercept=medianBC, lty=2)+
      theme(axis.text = element_text(size=11, color="black"),axis.title=element_text(size=15,color="black"))
    return(p)
}

#define a function to illustrate the basic distribution of BCcigars
BCcigars = function(umi_table, sample_name, color) {
    umi_table_selected = umi_table[,c("cellBC","intBC","CIGAR")]
    umi_table_selected$BCcigar = paste(umi_table_selected$intBC,umi_table_selected$CIGAR,sep=":")

    BCcigar = data.frame(table(umi_table_selected$BCcigar))
    BCcigar = BCcigar[order(BCcigar$Freq,decreasing = T),]
    sum_BCcigar = sum(BCcigar$Freq)
   
   BCcigar$sample = sample_name

    MajorBCcigar = filter(BCcigar, cumsum(Freq)<sum_BCcigar*0.95)
    MajorBCcigar = nrow(MajorBCcigar)
    shannonIndex = round(shannon(BCcigar$Freq),2)
    cat(c("The number of effective BCcigars is ",MajorBCcigar, ", The shannox index for BCcigar is", shannonIndex, "\n"))
    
    selected_BCcigar = BCcigar[1:40,]
    cigarIndel = as.character(selected_BCcigar$Var1)
    selected_BCcigar$Var1 = sapply(cigarIndel, function(eachIndel) { ifelse(nchar(eachIndel)>34,  substr(eachIndel,1, 34), eachIndel)}  )
    p = ggplot(selected_BCcigar, aes(x = factor(Var1,levels = unique(Var1)), y = 100*Freq/sum_BCcigar))+theme_classic()+
      scale_y_continuous(expand=c(0,0))+
      geom_col(fill = color)+theme(axis.text.x = element_text(angle = 90))+xlab("")+ylab("Percentage")+
       theme(axis.text = element_text(size=12, color="black"),axis.title=element_text(size=15,color="black"))
    result = list(BCcigar, p, MajorBCcigar, shannonIndex)
    return(result)
}


#need to use a for loop to handle multiple samples? How to save the intermediate results for further calling

#define a function to illustrate the basic distribution of CIGARs
cigar = function(umi_table, sample_name, color) {
    cigar = data.frame(table(umi_table$CIGAR))
    cigar = cigar[order(cigar$Freq,decreasing = T),]

    mut_rate = 1-cigar$Freq[cigar$Var1=="250M"]/sum(cigar$Freq) 
    mut_rate = round(mut_rate,2)
    cat(c("the total mutation rate for samples is ", mut_rate,". \n"))

    sum_cigar = sum(cigar$Freq)
    cigar$sample = sample_name

    Majorcigar = filter(cigar, cumsum(Freq)<sum_cigar*0.95)
    Majorcigar = nrow(Majorcigar)
    shannonIndex = round(shannon(cigar$Freq),2)
    cat(c("The number of effective cigars is ",Majorcigar," The shannox index for ", sample_name, " cigar is", shannonIndex, "\n"))
    selected_cigar = cigar[1:20,]
    cigarIndel = as.character(selected_cigar$Var1)
    selected_cigar$Var1 = sapply(cigarIndel, function(eachIndel) { ifelse(nchar(eachIndel)>20,  substr(eachIndel,1, 20), eachIndel)}  )
    p = ggplot( selected_cigar, aes(x = factor(Var1,levels = unique(Var1)), y = Freq*100/sum_cigar))+theme_classic()+
      scale_y_continuous(expand=c(0,0), limits=c(0,40))+ylab("Percentage")+
      geom_col(fill = color)+theme(axis.text.x = element_text(angle = 90))+xlab("")+
      ggtitle(sample_name)+
       theme(axis.text = element_text(size=12, color="black"),axis.title=element_text(size=15,color="black"))
    result = list(cigar, p, Majorcigar, shannonIndex, mut_rate)
    return(result)
}

#define a function to illustrate the basic distribution of intBCs
intBC = function(umi_table, sample_name) {
    intBC = data.frame(table(umi_table$intBC))
    intBC = intBC[order(intBC$Freq,decreasing = T),]
    intBC$sample = sample_name
    sum_intBC = sum(intBC$Freq)
    MajorintBC = filter(intBC, cumsum(Freq)<sum_intBC*0.95)
    MajorintBC = nrow( MajorintBC)
    shannonIndex = round(shannon(intBC$Freq),2)
    cat(c("The number of effective intBC is ", MajorintBC," The shannox index for intBC is", shannonIndex,".\n"))
    result = list(intBC, MajorintBC, shannonIndex)
    return(result)
}


#define function to handle CIGAR mutation
cigarnmo = function(umi_table, sample_name, color) {
    cigar = data.frame(table(umi_table$CIGAR))
    cigar = cigar[order(cigar$Freq,decreasing = T),]

    mut_rate = 1-cigar$Freq[cigar$Var1=="250M"]/sum(cigar$Freq) 
    mut_rate = round(mut_rate,2)
    cat(c("the total mutation rate for samples is ", mut_rate,". \n"))

    sum_cigar = sum(cigar$Freq)
    cigar$sample = sample_name

    Majorcigar = filter(cigar, cumsum(Freq)<sum_cigar*0.95)
    Majorcigar = nrow(Majorcigar)
    shannonIndex = round(shannon(cigar$Freq),2)
    cat(c("The number of effective cigars is ",Majorcigar," The shannox index for ", sample_name, " cigar is", shannonIndex, "\n"))
    selected_cigar = cigar[1:20,]
    cigarIndel = as.character(selected_cigar$Var1)
    selected_cigar$Var1 = sapply(cigarIndel, function(eachIndel) { ifelse(nchar(eachIndel)>20,  substr(eachIndel,1, 20), eachIndel)}  )
    p = ggplot( selected_cigar, aes(x = factor(Var1,levels = unique(Var1)), y = Freq*100/sum_cigar))+theme_classic()+
      scale_y_continuous(expand=c(0,0), limits=c(0,40))+ylab("Relative frequency (%)")+
      geom_col(fill = color)+theme(axis.text.x = element_text(angle = 90))+xlab("")+
       theme(axis.text = element_text(size=12, color="black"),axis.title=element_text(size=15,color="black"))
    result = list(cigar, p, Majorcigar, shannonIndex, mut_rate)
    return(result)
}


#this function can only apply to individual target sites
site_indel = function(site, color, ylimit=16) {
        site = gsub("[ATCG]", "",  site) %>% str_replace_all(fixed("["), "") %>% str_replace_all(fixed("]"), "")    #remove unecessary strings
        site = site[site!=""]     #remove those reads that could not call the indels
        #r1 = sapply(r1, function(x) ifelse(length(str_split(x, fixed("["))[[1]])>1, str_split(x, fixed("["))[[1]][2], str_split(x, fixed("["))[[1]]))   #running too slow, remove those indels call twice in different targets
        site = data.frame(table(site))
        site = site[order(site$Freq,decreasing = T),]
        sum_indel = sum(site$Freq)
        site$freq = site$Freq/sum_indel
        mutation_rate = round(1-site$freq[1], 3)  #this mutation rate is higher since some un-called indels are also incorporated
        diversity_index = round(shannon(site$Freq), 2) #index for both mutated alleles and non-mutated alleles
        cat(c("the mutation rate is ", mutation_rate,", and the shannon index is ", diversity_index, ". " )) #consider the true shannon index
    
        site = site[-1,]  #only consider those reads show mutations
        sum_indel = sum(site$Freq) #calculate the sum of indels again, since the most frequent one have been removed
        Major_indel = nrow(filter(site, cumsum(Freq)<sum_indel*0.95))  #effective indels, only consider muatated indels
        
        cat(c("Mutated alleles only. The number of effective indels is ",Major_indel, "\n"))

        site$site = factor(site$site,levels = site$site)
        options(repr.plot.width= 8, repr.plot.height=8)
        site_top20 = head(site, 20)  #be careful, indel frequency still consider all alleles, including the non-mutated alleles
        site_top20$freq = 100*site_top20$freq
        p = ggplot(site_top20, aes(site,freq))+geom_bar(stat="identity", fill = color )+theme_classic()+scale_y_continuous(expand=c(0,0), limits = c(0, ylimit))+
            theme(axis.text.x = element_text(angle = 90))+ggtitle(paste0(Major_indel, " unique indels"))+xlab("")+ylab("Relative indel frquency (%)")+
             theme(axis.text = element_text(size=11, color="black"),axis.title=element_text(size=15,color="black"))
        return(p)
    }

#using a percent threshold to remove possible intBCs from sequencing errors
filterIntBC = function(sample, filter_percent) {
    #filename1 = paste0('../targetProcessData/', sample, '_preprocess_pipeline/umi_table.csv')
    umi_table  = sample
    intBC = data.frame(table(umi_table$intBC[nchar(umi_table$intBC)==14]))
    intBC = intBC[order(intBC$Freq, decreasing=T),]
    intBC1 = intBC[cumsum(intBC$Freq)<filter_percent*sum(intBC$Freq),]
    colnames(intBC1)=c("intBC","Freq")
    candidate_barcodes = intBC1$intBC
    umi_table_filter = umi_table %>% filter(intBC %in% candidate_barcodes) 
    #filename2 = paste0('/data1/home/gdpeng/chengchen/dualproject/targetProcessData/', sample, '_preprocess_pipeline/umi_table_filtered.csv')
    #write.csv(umi_table_filter, filename2)
    return(umi_table_filter)
    }


#===============================================================================
# For mutation plots similar to DARLIN and CARLIN
#input: character, a single CIGAR value
#output: data.frame, three columns, mutation start site, mutation end site, mutation type
parse_cigar = function(cigar) 
{
    operations <- strsplit(gsub("(\\D)", "\\1 ", cigar), " ")[[1]]
    position <- 0
    results <- list()
    for (operation in operations) {
        length <- as.numeric(gsub("\\D", "", operation))
        type <- gsub("\\d", "", operation)
        if (type %in% c("D", "I")) {
            results[[length(results) + 1]] <- list(start = position + 
                1, end = position + length, type = type)
        }
        if (type %in% c("M", "D")) {
            position <- position + length
        }
    }
    
    results = do.call(rbind, lapply(results, data.frame, stringsAsFactors=FALSE))
    return(results)
}


#considering using UMI fraction instead of cell fraction, cell fraction are harder to evaluate, 
# and not represent the true mutation sceniario
#considering that a cell may have two intBCs, each intBC harbour the same mutation or different mutation
#it will be extremely difficult for further distinction different mutations from different cells.

#input: umi_table data.frame, read.csv without specifying function row.names=1
#input: top numeric, most common indels, top=0 means no selection
#input: df_ref dataframe, ref-sequence annotation
#output: data.frame, with 7 columns: start	end	type	n_umi	fraction	log10fraction	linewidth
umi_table_to_df_mut = function(umi_table, df_ref, top=200) {
    nUMI = nrow(umi_table)
    cigar_umi = umi_table %>% group_by(CIGAR) %>% summarise(n_umi = n_distinct(X)) %>% arrange(desc(n_umi)) 
    # Test the function
    parse_results = data.frame(start=double(), end=double(), type=character(), n_umi=double())
    for (i in 1:nrow(cigar_umi)) {
        cigar = cigar_umi$CIGAR[i]
        if (grepl('[ID]', cigar)) {
            result = parse_cigar(cigar)

            result$n_umi = cigar_umi$n_umi[i]
            parse_results = rbind(parse_results, result)
            }
    }

    cigar_pattern = parse_results %>% group_by(start, end, type) %>% reframe(n_umi = sum(n_umi)) %>% 
       arrange(desc(n_umi)) %>% mutate(fraction = n_umi/nUMI) 

    #only check the top 200 mutations or all mutations?    
    df_mut = cigar_pattern %>% filter(start>=df_ref$start[1] & end<df_ref$end[nrow(df_ref)]) 
    if (top > 0) {
        df_mut = head(df_mut, top) }
    df_mut$log10fraction = log10(df_mut$fraction)
    df_mut$linewidth = -3/(log10(df_mut$fraction)) #3 was selected randomly accordint to the linewidth
    
    #some linewidth is extremely large, single only a few mutations occurs in that intBCs, need to constrain its linewidth
    df_mut$linewidth = ifelse(df_mut$linewidth>10, 10, df_mut$linewidth) #nearly 50%
    
    # df_mut = df_mut %>% filter(start>=df_ref$start[1]) #some mutation are outside of flanking region
    return(df_mut)
}


#create only one sector and then distinguish different fragment using circos.rect(), circos.arrow, and circos.line

# Assuming you have a data frame `df_ref` with columns 'region', 'start', 'end'
# And a data frame `df_mut` with columns 'type', 'start', 'end' n_umi	fraction	log10fraction	linewidth

library(circlize)
library(ComplexHeatmap)
mutation_circos_plot = function(df_ref, df_mut) {
        
    total_rects = nrow(df_ref)

    # Define the genomic information for the circos plot
    circos.clear()
    circos.par(start.degree = 90)
    # circos.par(gap.degree = 0.5, cell.padding = c(0.02, 0.02, 0.02, 0.02))

    circos.initialize('backbone', xlim = c(df_ref$start[1],df_ref$end[total_rects]))

    circos.track(ylim = c(0, 1), 
                 # bg.col = col,
        bg.col = c('white'), 
        bg.border = NA, track.height = 0.2)

    #plot the regions
    for (i in 1: (total_rects-1)) {
        circos.rect(xleft=df_ref$start[i], ybottom=0, xright=df_ref$end[i], ytop=1,  col=colors[df_ref$region[i]])
    }

    #show the direction of fragments
    circos.arrow(
        x1=df_ref$start[total_rects],
        x2=df_ref$end[total_rects],
        y=get.cell.meta.data("ycenter"),
        width = get.cell.meta.data("yrange"),
        arrow.head.length = mm_x(10),
        arrow.head.width = get.cell.meta.data("yrange")*1.2,
        col='#808080'
        )

    #add labels
    for (i in 1: total_rects) {
        circos.text(x=(df_ref$start[i]+df_ref$end[i])/2, y=0.5, df_ref$region[i]) #adjust the font size
    }

    # Add the links for each mutation
    for (i in 1:nrow(df_mut)) {
        circos.link('backbone',
                    df_mut$start[i], 
                    'backbone', 
                    df_mut$end[i], 
                    col = ifelse(df_mut$type[i] == "D", "#fe0000", "#0000fe"), #deletion or insertions
                    lwd = df_mut$linewidth[i]) #linewidth according to frequency

    }

    #add legends
    draw(Legend(at = c('Deletion', 'Insertion'), type='points', legend_gp=gpar(col=c("#fe0000", "#0000fe")),
        title_position = "topleft", title = "Links"), 
         # x = unit(1, "npc") - unit(3, "mm"), y = unit(8, "mm"),
         x = unit(0.95, "npc"), y = unit(0.05, "npc"),
         just=c('right', 'bottom'))

    draw(Legend(at=(-2), type='lines', legend_gp=gpar(col='#81007f', lwd=1.5), title_position='topleft', title='log10 UMI\n Fraction'), 
        # x = unit(6, "mm"), y = unit(8, "mm"),
         x = unit(0.05, "npc"), y = unit(0.05, "npc"),
        just=c('left', 'bottom'))

    # add a title for the plot to specify the number of mutations
    n = nrow(df_mut)
    title(paste0('Total number of mutation types is ', n))
    circos.clear() }


# #generate UMI barcode fraction for different types of mutation
# #deal with each target site individually


process_row_r1 <- function(i, umi_table, type) {
  r1 = unlist(umi_table$r1[i])
  r2 = unlist(umi_table$r2[i])
  r3 = unlist(umi_table$r3[i])
  r4 = unlist(umi_table$r4[i])
  
  if (length(r1)==1 && r1[1]=='None') {
    return('Unedited')
  } else if (length(r1)==1 && any(grepl('D', r1)) && !(r1[1] %in% c(r2,r3,r4))) {
    return('Intrasite_deletion') #whether r1 mutation event is in r2/r3/r4
  } else if (length(r1)==1 && any(grepl('D', r1)) && (r1[1] %in% c(r2,r3,r4))) {
    return('Intersite_deletion')
  } else if (length(r1)>1 && (r1[length(r1)] %in% c(r2,r3,r4))) {
    return('Intersite_indel')
  } else if (length(r1)>1 && !(r1[length(r1)] %in% c(r2,r3,r4))) {
    return('Intrasite_indel')
  } else if (length(r1)==1 && any(grepl('I', r1))) {
    return('Insertion')
  } else {
    return(umi_table$r1_mut_type[i])
  }
}

process_row_r2 <- function(i, umi_table, type) {
  r1 = unlist(umi_table$r1[i])
  r2 = unlist(umi_table$r2[i])
  r3 = unlist(umi_table$r3[i])
  r4 = unlist(umi_table$r4[i])

 if (length(r2)==1 && r2[1]=='None') {
    return('Unedited')
  } else if (length(r2)==1 && any(grepl('D', r2)) && !(r2[1] %in% c(r1,r3,r4))) {
    return('Intrasite_deletion') 
  } else if (length(r2)==1 && any(grepl('D', r2)) && (r2[1] %in% c(r1,r3,r4))) {
    return('Intersite_deletion')
  } else if (length(r2)>1 && (r2[length(r2)] %in% c(r1,r3,r4))) {
    return('Intersite_indel')
  } else if (length(r2)>1 && !(r2[length(r2)] %in% c(r1,r3,r4))) {
    return('Intrasite_indel')
  } else if (length(r2)==1 && any(grepl('I', r2))) {
    return('Insertion')
  } else {
    return(umi_table$r2_mut_type[i])
  }
}

process_row_r3 <- function(i, umi_table, type) {
  r1 = unlist(umi_table$r1[i])
  r2 = unlist(umi_table$r2[i])
  r3 = unlist(umi_table$r3[i])
  r4 = unlist(umi_table$r4[i])
  
  if (length(r3)==1 && r3[1]=='None') {
    return('Unedited')
  } else if (length(r3)==1 && any(grepl('D', r3)) && !(r3[1] %in% c(r1,r2,r4))) {
    return('Intrasite_deletion')
  } else if (length(r3)==1 && any(grepl('D', r3)) && (r3[1] %in% c(r1,r2,r4))) {
    return('Intersite_deletion')
  } else if (length(r3)>1 && (r3[length(r3)] %in% c(r1,r2,r4))) {
    return('Intersite_indel')
  } else if (length(r3)>1 && !(r3[length(r3)] %in% c(r1,r2,r4))) {
    return('Intrasite_indel')
  } else if (length(r3)==1 && any(grepl('I', r3))) {
    return('Insertion')
  } else {
    return(umi_table$r3_mut_type[i])
  }
}

process_row_r4 <- function(i, umi_table, type) {
  r1 = unlist(umi_table$r1[i])
  r2 = unlist(umi_table$r2[i])
  r3 = unlist(umi_table$r3[i])
  r4 = unlist(umi_table$r4[i])
  
  if (length(r4)==1 && r4[1]=='None') {
    return('Unedited')
  } else if (length(r4)==1 && any(grepl('D', r4)) && !(r4[1] %in% c(r1, r2,r3))) {
    return('Intrasite_deletion')
  } else if (length(r4)==1 && any(grepl('D', r4)) && (r4[1] %in% c(r1,r2,r3))) {
    return('Intersite_deletion')
  } else if (length(r4)>1 && (r4[length(r4)] %in% c(r1,r2,r3))) {
    return('Intersite_indel')
  } else if (length(r4)>1 && !(r4[length(r4)] %in% c(r1,r2,r3))) {
    return('Intrasite_indel')
  } else if (length(r4)==1 && any(grepl('I', r4))) {
    return('Insertion')
  } else {
    return(umi_table$r4_mut_type[i])
  }
}



# process_row_r1 <- function(i, umi_table, type) {
#   r1 = unlist(umi_table$r1[i])
#   r2 = unlist(umi_table$r2[i])
#   r3 = unlist(umi_table$r3[i])
#   r4 = unlist(umi_table$r4[i])
  
#   if (type == 'Unedited' && length(r1)==1 && r1[1]=='None') {
#     return('Unedited')
#   } else if (type == 'Intrasite_deletion' && length(r1)==1 && any(grepl('D', r1)) && !(r1[1] %in% c(r2,r3,r4))) {
#     return('Intrasite_deletion') #whether r1 mutation event is in r2/r3/r4
#   } else if (type == 'Intersite_deletion' && length(r1)==1 && any(grepl('D', r1)) && (r1[1] %in% c(r2,r3,r4))) {
#     return('Intersite_deletion')
#   } else if (type == 'Intersite_indel' && length(r1)>1 && (r1[length(r1)] %in% c(r2,r3,r4))) {
#     return('Intersite_indel')
#   } else if (type == 'Intrasite_indel' && length(r1)>1 && !(r1[length(r1)] %in% c(r2,r3,r4))) {
#     return('Intrasite_indel')
#   } else if (type == 'Insertion' && length(r1)==1 && any(grepl('I', r1))) {
#     return('Insertion')
#   } else {
#     return(umi_table$r1_mut_type[i])
#   }
# }

# process_row_r2 <- function(i, umi_table, type) {
#   r1 = unlist(umi_table$r1[i])
#   r2 = unlist(umi_table$r2[i])
#   r3 = unlist(umi_table$r3[i])
#   r4 = unlist(umi_table$r4[i])

#  if (type == 'Unedited' && length(r2)==1 && r2[1]=='None') {
#     return('Unedited')
#   } else if (type == 'Intrasite_deletion' && length(r2)==1 && any(grepl('D', r2)) && !(r2[1] %in% c(r1,r3,r4))) {
#     return('Intrasite_deletion') 
#   } else if (type == 'Intersite_deletion' && length(r2)==1 && any(grepl('D', r2)) && (r2[1] %in% c(r1,r3,r4))) {
#     return('Intersite_deletion')
#   } else if (type == 'Intersite_indel' && length(r2)>1 && (r2[length(r2)] %in% c(r1,r3,r4))) {
#     return('Intersite_indel')
#   } else if (type == 'Intrasite_indel' && length(r2)>1 && !(r2[length(r2)] %in% c(r1,r3,r4))) {
#     return('Intrasite_indel')
#   } else if (type == 'Insertion' && length(r2)==1 && any(grepl('I', r2))) {
#     return('Insertion')
#   } else {
#     return(umi_table$r2_mut_type[i])
#   }
# }

# process_row_r3 <- function(i, umi_table, type) {
#   r1 = unlist(umi_table$r1[i])
#   r2 = unlist(umi_table$r2[i])
#   r3 = unlist(umi_table$r3[i])
#   r4 = unlist(umi_table$r4[i])
  
#   if (type == 'Unedited' && length(r3)==1 && r3[1]=='None') {
#     return('Unedited')
#   } else if (type == 'Intrasite_deletion' && length(r3)==1 && any(grepl('D', r3)) && !(r3[1] %in% c(r1,r2,r4))) {
#     return('Intrasite_deletion')
#   } else if (type == 'Intersite_deletion' && length(r3)==1 && any(grepl('D', r3)) && (r3[1] %in% c(r1,r2,r4))) {
#     return('Intersite_deletion')
#   } else if (type == 'Intersite_indel' && length(r3)>1 && (r3[length(r3)] %in% c(r1,r2,r4))) {
#     return('Intersite_indel')
#   } else if (type == 'Intrasite_indel' && length(r3)>1 && !(r3[length(r3)] %in% c(r1,r2,r4))) {
#     return('Intrasite_indel')
#   } else if (type == 'Insertion' && length(r3)==1 && any(grepl('I', r3))) {
#     return('Insertion')
#   } else {
#     return(umi_table$r3_mut_type[i])
#   }
# }

# process_row_r4 <- function(i, umi_table, type) {
#   r1 = unlist(umi_table$r1[i])
#   r2 = unlist(umi_table$r2[i])
#   r3 = unlist(umi_table$r3[i])
#   r4 = unlist(umi_table$r4[i])
  
#   if (type == 'Unedited' && length(r4)==1 && r4[1]=='None') {
#     return('Unedited')
#   } else if (type == 'Intrasite_deletion' && length(r4)==1 && any(grepl('D', r4)) && !(r4[1] %in% c(r1, r2,r3))) {
#     return('Intrasite_deletion')
#   } else if (type == 'Intersite_deletion' && length(r4)==1 && any(grepl('D', r4)) && (r4[1] %in% c(r1,r2,r3))) {
#     return('Intersite_deletion')
#   } else if (type == 'Intersite_indel' && length(r4)>1 && (r4[length(r4)] %in% c(r1,r2,r3))) {
#     return('Intersite_indel')
#   } else if (type == 'Intrasite_indel' && length(r4)>1 && !(r4[length(r4)] %in% c(r1,r2,r3))) {
#     return('Intrasite_indel')
#   } else if (type == 'Insertion' && length(r4)==1 && any(grepl('I', r4))) {
#     return('Insertion')
#   } else {
#     return(umi_table$r4_mut_type[i])
#   }
# }



#DARLIN/CARLIN's way to define the mutation type
# for each target site, consider the last mutation events, and perform pairwise comparison to 
# check whether the early mutation events also appeared at target site afterwards
#if not [none] and 'NA', split the target site according to the number of brackets
# input: umi_table, a data.frame generated from Cassiopeia preprocessing. Note: This umi_table store one UMI per row
# output: a umi_fraction, a data.frame stores the umi Freq for each type of target site

mutation_umi_fraction = function(umi_table) {

    # preprocessing
    umi_table = umi_table %>% select(X, r1, r2, r3, r4)
    umi_table$r1 = gsub('[ATCG]', '', umi_table$r1)
    umi_table$r2 = gsub('[ATCG]', '', umi_table$r2)
    umi_table$r3 = gsub('[ATCG]', '', umi_table$r3)
    umi_table$r4 = gsub('[ATCG]', '', umi_table$r4)

    umi_table$r1 = ifelse(umi_table$r1=='', 'unknown', umi_table$r1)
    umi_table$r2 = ifelse(umi_table$r2=='', 'unknown', umi_table$r2)
    umi_table$r3 = ifelse(umi_table$r3=='', 'unknown', umi_table$r3)
    umi_table$r4 = ifelse(umi_table$r4=='', 'unknown', umi_table$r4)

    #split mutation events within a target
    umi_table$r1 = sub('^\\[', '', umi_table$r1)
    umi_table$r2 = sub('^\\[', '', umi_table$r2)
    umi_table$r3 = sub('^\\[', '', umi_table$r3)
    umi_table$r4 = sub('^\\[', '', umi_table$r4)

    umi_table$r1 = sub('\\]$', '', umi_table$r1)
    umi_table$r2 = sub('\\]$', '', umi_table$r2)
    umi_table$r3 = sub('\\]$', '', umi_table$r3)
    umi_table$r4 = sub('\\]$', '', umi_table$r4)

    umi_table$r1 = strsplit(umi_table$r1, "\\]\\[")
    umi_table$r2 = strsplit(umi_table$r2, "\\]\\[")
    umi_table$r3 = strsplit(umi_table$r3, "\\]\\[")
    umi_table$r4 = strsplit(umi_table$r4, "\\]\\[")

    # assign mutation event information for each target
    
    umi_table$r1_mut_type = 'unknown'
    umi_table$r2_mut_type = 'unknown'
    umi_table$r3_mut_type = 'unknown'
    umi_table$r4_mut_type = 'unknown'
    
    umi_table$r1_mut_type = sapply(1:nrow(umi_table), process_row_r1, umi_table=umi_table)
    umi_table$r2_mut_type = sapply(1:nrow(umi_table), process_row_r2, umi_table=umi_table)
    umi_table$r3_mut_type = sapply(1:nrow(umi_table), process_row_r3, umi_table=umi_table)
    umi_table$r4_mut_type = sapply(1:nrow(umi_table), process_row_r4, umi_table=umi_table)    


    # Jan 7th, 2024
    r1_umi_fraction = table(umi_table$r1_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
      filter(type!='unknown')  %>% mutate(target = 'T1')
    r2_umi_fraction = table(umi_table$r2_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
      filter(type!='unknown')  %>% mutate(target = 'T2')
    r3_umi_fraction = table(umi_table$r3_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
      filter(type!='unknown')  %>% mutate(target = 'T3')
    r4_umi_fraction = table(umi_table$r4_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
      filter(type!='unknown')  %>% mutate(target = 'T4')
    umi_fraction = rbind(r1_umi_fraction, r2_umi_fraction, r3_umi_fraction, r4_umi_fraction)
    umi_fraction$type = factor(umi_fraction$type, levels = c('Unedited', 'Intrasite_deletion', 'Intersite_deletion', 
    'Insertion','Intrasite_indel', 'Intersite_indel'))
    return(umi_fraction)
}


# mutation_umi_fraction = function(umi_table) {
    
#     # preprocessing
#     umi_table = umi_table %>% select(X, r1, r2, r3, r4)
#     umi_table$r1 = gsub('[ATCG]', '', umi_table$r1)
#     umi_table$r2 = gsub('[ATCG]', '', umi_table$r2)
#     umi_table$r3 = gsub('[ATCG]', '', umi_table$r3)
#     umi_table$r4 = gsub('[ATCG]', '', umi_table$r4)

#     umi_table$r1 = ifelse(umi_table$r1=='', 'unknown', umi_table$r1)
#     umi_table$r2 = ifelse(umi_table$r2=='', 'unknown', umi_table$r2)
#     umi_table$r3 = ifelse(umi_table$r3=='', 'unknown', umi_table$r3)
#     umi_table$r4 = ifelse(umi_table$r4=='', 'unknown', umi_table$r4)
    
#     #split mutation events within a target
#     umi_table$r1 = sub('^\\[', '', umi_table$r1)
#     umi_table$r2 = sub('^\\[', '', umi_table$r2)
#     umi_table$r3 = sub('^\\[', '', umi_table$r3)
#     umi_table$r4 = sub('^\\[', '', umi_table$r4)

#     umi_table$r1 = sub('\\]$', '', umi_table$r1)
#     umi_table$r2 = sub('\\]$', '', umi_table$r2)
#     umi_table$r3 = sub('\\]$', '', umi_table$r3)
#     umi_table$r4 = sub('\\]$', '', umi_table$r4)

#     umi_table$r1 = strsplit(umi_table$r1, "\\]\\[")
#     umi_table$r2 = strsplit(umi_table$r2, "\\]\\[")
#     umi_table$r3 = strsplit(umi_table$r3, "\\]\\[")
#     umi_table$r4 = strsplit(umi_table$r4, "\\]\\[")
    
#     # assign mutation event information for each target
#     types = c('Unedited', 'Intrasite_deletion', 'Intersite_deletion', 'Intersite_indel', 'Intrasite_indel', 'Insertion')

#     umi_table$r1_mut_type = 'unknown'
#     umi_table$r2_mut_type = 'unknown'
#     umi_table$r3_mut_type = 'unknown'
#     umi_table$r4_mut_type = 'unknown'
#     for (type in types) {
#         umi_table$r1_mut_type = sapply(1:nrow(umi_table), process_row_r1, umi_table=umi_table, type=type)
#         umi_table$r2_mut_type = sapply(1:nrow(umi_table), process_row_r2, umi_table=umi_table, type=type)
#         umi_table$r3_mut_type = sapply(1:nrow(umi_table), process_row_r3, umi_table=umi_table, type=type)
#         umi_table$r4_mut_type = sapply(1:nrow(umi_table), process_row_r4, umi_table=umi_table, type=type)    
#     }

#     # Jan 7th, 2024
#     r1_umi_fraction = table(umi_table$r1_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
#       filter(type!='unknown')  %>% mutate(target = 'T1')
#     r2_umi_fraction = table(umi_table$r2_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
#       filter(type!='unknown')  %>% mutate(target = 'T2')
#     r3_umi_fraction = table(umi_table$r3_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
#       filter(type!='unknown')  %>% mutate(target = 'T3')
#     r4_umi_fraction = table(umi_table$r4_mut_type) %>% as.data.frame() %>% rename(type=Var1, Freq=Freq) %>% 
#       filter(type!='unknown')  %>% mutate(target = 'T4')
#     umi_fraction = rbind(r1_umi_fraction, r2_umi_fraction, r3_umi_fraction, r4_umi_fraction)
#     umi_fraction$type = factor(umi_fraction$type, levels = c('Unedited', 'Intrasite_deletion', 'Intersite_deletion', 
#     'Insertion','Intrasite_indel', 'Intersite_indel'))
#     return(umi_fraction)
# }


#try to compare strings to evaluate the deletion status
#input: umi_table, dataframe, generated from cassiopeia process
#output: inter_deletion_fraction, dataframe
intersite_deletion_fraction = function(umi_table) {
    umi_table = umi_table %>% select(X, r1, r2, r3, r4)
    umi_table$r1 = gsub('[ATCG]', '', umi_table$r1)
    umi_table$r2 = gsub('[ATCG]', '', umi_table$r2)
    umi_table$r3 = gsub('[ATCG]', '', umi_table$r3)
    umi_table$r4 = gsub('[ATCG]', '', umi_table$r4)

    umi_table$r1_r2del =  sapply(1:nrow(umi_table), function(i) ifelse(umi_table$r1[i] != "[None]" & umi_table$r1[i] !="", grepl(umi_table$r1[i], umi_table$r2[i], fixed=TRUE), FALSE))
    umi_table$r1_r3del =  sapply(1:nrow(umi_table), function(i) ifelse(umi_table$r1[i] != "[None]" & umi_table$r1[i] !="", grepl(umi_table$r1[i], umi_table$r3[i], fixed=TRUE), FALSE))
    umi_table$r1_r4del =  sapply(1:nrow(umi_table), function(i) ifelse(umi_table$r1[i] != "[None]" & umi_table$r1[i] !="", grepl(umi_table$r1[i], umi_table$r4[i], fixed=TRUE), FALSE))
    umi_table$r2_r3del =  sapply(1:nrow(umi_table), function(i) ifelse(umi_table$r2[i] != "[None]" & umi_table$r2[i] !="", grepl(umi_table$r2[i], umi_table$r3[i], fixed=TRUE), FALSE))
    umi_table$r2_r4del =  sapply(1:nrow(umi_table), function(i) ifelse(umi_table$r2[i] != "[None]" & umi_table$r2[i] !="", grepl(umi_table$r2[i], umi_table$r4[i], fixed=TRUE), FALSE)) 
    umi_table$r3_r4del =  sapply(1:nrow(umi_table), function(i) ifelse(umi_table$r3[i] != "[None]" & umi_table$r3[i] !="", grepl(umi_table$r3[i], umi_table$r4[i], fixed=TRUE), FALSE)) 

    #it's already known there is no r1-r4 deletion
    #r1_r3 deletion will also cover r1_r2 deletion and r2_r3 deletion, try to correct it
    umi_table$r1_r3del = ifelse(umi_table$r1_r4del & umi_table$r1_r3del, FALSE, umi_table$r1_r3del)
    umi_table$r1_r2del = ifelse(umi_table$r1_r4del & umi_table$r1_r2del, FALSE, umi_table$r1_r2del)  #r1_r3 del overlap with r1_r2del and r2_r3del
    umi_table$r2_r3del = ifelse(umi_table$r1_r4del & umi_table$r2_r3del, FALSE, umi_table$r2_r3del)
    umi_table$r2_r4del = ifelse(umi_table$r1_r4del & umi_table$r2_r4del, FALSE, umi_table$r2_r4del)                             
    umi_table$r2_r3del = ifelse(umi_table$r1_r4del & umi_table$r2_r3del, FALSE, umi_table$r2_r3del)
    umi_table$r3_r4del = ifelse(umi_table$r1_r4del & umi_table$r3_r4del, FALSE, umi_table$r3_r4del) 

    umi_table$r1_r2del = ifelse(umi_table$r1_r3del & umi_table$r1_r2del, FALSE, umi_table$r1_r2del)  #r1_r3 del overlap with r1_r2del and r2_r3del
    umi_table$r2_r3del = ifelse(umi_table$r1_r3del & umi_table$r2_r3del, FALSE, umi_table$r2_r3del)

    umi_table$r2_r3del = ifelse(umi_table$r2_r4del & umi_table$r2_r3del, FALSE, umi_table$r2_r3del)
    umi_table$r3_r4del = ifelse(umi_table$r2_r4del & umi_table$r3_r4del, FALSE, umi_table$r3_r4del) 

    #allele_level calculation of mutation levels, 
    #however, here, this calculation did not consider the information of readCounts, 
    #as UMI are more accurate estimation of actual molecular diversity, which is difficult to achieve during bulk sequencing.
    r1_r2ratio = sum(umi_table$r1_r2del)/nrow(umi_table)
    r2_r3ratio = sum(umi_table$r2_r3del)/nrow(umi_table) #actually this function only consider the mutation fraction in the allele_level
    r3_r4ratio = sum(umi_table$r3_r4del)/nrow(umi_table) #Not at the UMI level
    r1_r3ratio = sum(umi_table$r1_r3del)/nrow(umi_table)
    r2_r4ratio = sum(umi_table$r2_r4del)/nrow(umi_table)
    r1_r4ratio = sum(umi_table$r1_r4del)/nrow(umi_table)

    r1ratio = sum(grepl("[DI]", umi_table$r1))/nrow(umi_table) - r1_r2ratio - r1_r3ratio - r1_r4ratio
    r2ratio = sum(grepl("[DI]", umi_table$r2))/nrow(umi_table) - r1_r2ratio - r2_r3ratio - r1_r3ratio - r2_r4ratio - r1_r4ratio
    r3ratio = sum(grepl("[DI]", umi_table$r3))/nrow(umi_table) - r2_r3ratio - r3_r4ratio - r1_r3ratio - r2_r4ratio - r1_r4ratio
    r4ratio = sum(grepl("[DI]", umi_table$r4))/nrow(umi_table) - r3_r4ratio - r2_r4ratio -r1_r4ratio 

    #r2_r4ratio is too low to show
    ratio = c(r1ratio, r2ratio, r3ratio, r4ratio, r1_r2ratio, r2_r3ratio, r3_r4ratio, r1_r3ratio, r2_r4ratio, r1_r4ratio)*100
    state = c("T1", "T2", "T3", "T4", "T1-2", "T2-3", "T3-4", "T1-3", "T2-4", "T1-4")
    state = factor(state, levels=state)
    type = c("Intrasite", "Intrasite", "Intrasite", "Intrasite", "Inter-2-sites", "Inter-2-sites", "Inter-2-sites", "Inter-3-sites", "Inter-3-sites", "Inter-4-sites")
    type = factor(type, levels = c("Intrasite", "Inter-2-sites", "Inter-3-sites", "Inter-4-sites"))
    inter_deletion_fraction = data.frame(state, ratio, type)
    return(inter_deletion_fraction)
}


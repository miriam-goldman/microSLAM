
validate_genes_input<-function(opt){
  s_id=opt$species_id
  verbose=opt$verbose
  put(paste("species id for out files is:",s_id),console = verbose)
    
  ##read in gene data
  if(isTRUE(!is.na(opt$genes_copynumber_file))){
    if(file_test("-f",opt$genes_copynumber_file)){
      gcopynumber<-fread(opt$genes_copynumber_file)
    }else{
      put("copy number file doesnt exist",console = verbose)
      stop()
    }
  }else{
    put("copy number file doesnt exist",console = verbose)
    stop()
  }
  
  if(isTRUE(!is.na(opt$genes_copynumber_file))){
    if(file_test("-f",opt$genes_depth_file)){
      gdepth<-fread(opt$genes_depth_file)
    }else{
      put("genes depth file doesnt exist",console = verbose)
    }
  }else{
    put("genes depth file doesnt exist",console = verbose)
  }
  
     put("gene files read in",console = verbose)
  
  
  if(isTRUE(!is.na(opt$genes_summary))){
    if(file_test("-f",opt$genes_summary)){
      if(verbose){
        put("reading in genes summary",console = verbose)
      }
      genes_summary<-fread(opt$genes_summary) %>% filter(marker_coverage > 0)
      genes_summary_used=TRUE
    }else{
      genes_summary<-NULL
      genes_summary_used=FALSE
      put("genes depth file doesnt exist",console = verbose)
    }
  }else{
    genes_summary<-NULL
    genes_summary_used=FALSE
    put("genes summary file doesnt exist",console = verbose)
  }
  
  if(isTRUE(!is.na(opt$centroid_prevalence_file))){
    if(file_test("-f",opt$centroid_prevalence_file)){
      centroid_prevalence_file<-fread(opt$centroid_prevalence_file)
      pangenome_used=TRUE
      centroid_prevalence_cutoff=opt$centroid_prevalence_cutoff
      colnames(centroid_prevalence_file)<-c("rep_gene_id","centroid_prevalence",	"centroid_ength")
      put(paste("pangenome used and read in centriod cutoff at",centroid_prevalence_cutoff),console = verbose)
    }else{
      centroid_prevalence_cutoff<-NA
      centroid_prevalence_file<-NULL
      pangenome_used=FALSE
    }
  }else{
    centroid_prevalence_cutoff<-NA
    centroid_prevalence_file<-NULL
    pangenome_used=FALSE
  }
     
     if(isTRUE(!is.na(opt$genes_info))){
       if(file_test("-f",opt$genes_info)){
         genes_info<-fread(opt$genes_info)
         genes_info<-genes_info %>% group_by(centroid_80) %>% summarize(n=n()) %>% mutate(singleton=ifelse(n==1,TRUE,FALSE))
         genes_info<-genes_info %>% right_join(centroid_prevalence_file,by=c("centroid_80"="rep_gene_id"))       
         }
       }
     
     
  if(isTRUE(!is.na(opt$metadata))){
    if(file_test("-f",opt$metadata)){
      metadata <- fread(opt$metadata) 
      metadata_overlaps<-any(metadata$sample_name %in% colnames(gcopynumber))
      if(!metadata_overlaps){
        put("samples in metadata do not match MIDAS output please make sure you metadata has sample names
            in a colmun labeled sample_name and binary phenotypes in a column labeled disease_status",console = verbose)
      }else{
        metadata_used=TRUE
        min_control_case_ratio=opt$min_control_case_ratio
        put(paste("using metadata to filter, filtering controls to case ratio at",min_control_case_ratio,"per gene"),console = verbose)
      }
    }else{
      min_control_case_ratio=opt$min_control_case_ratio
      metadata=NULL
      metadata_used=FALSE
    }
  }else{
    min_control_case_ratio=opt$min_control_case_ratio
    metadata=NULL
    metadata_used=FALSE
    if(verbose){
      put("not using metadata",console = verbose)
    }
    
  }
  if(isTRUE(!is.na(opt$GRM))){
    if(file_test("-f",opt$GRM)){
      GRM <- fread(opt$GRM,sep="\t",header=FALSE) 
      colnames(GRM)<-c("sample_name",GRM$V1)
      setindexv(GRM,'sample_name')
      GRM_overlaps<-any(colnames(GRM) %in% colnames(gcopynumber))
      if(!GRM_overlaps){
        put("samples in GRM do not match MIDAS output please make sure your GRM has sample names as the colmun names",console = verbose)
        stop()
      }else{
        GRM_used=TRUE
          put(paste("using GRM to filter sample list from GRM will be used"),console = verbose)
      }
    }
  }else{
    GRM=NULL
    GRM_used=FALSE
    put("not using GRM to filter",console = verbose)
    
  }
  samples_per_copynumber<-opt$number_of_samples_for_copynumber
  if(!is.numeric(samples_per_copynumber)){
    put("samples per copynumber invlaid use number or 0",console = verbose)
  }
  depth_cutoff<-opt$depth_cutoff
  if(!is.numeric(depth_cutoff)){
    put("depth cutoff invlaid use number or 0",console = verbose)
  }
  if(is.logical(opt$log_scale)){
    log_scale=opt$log_scale
  }else{
    log_scale=FALSE
  }
  if(is.logical(opt$mean_center)){
    mean_center=opt$mean_center
  }else{
    mean_center=FALSE
  }
  if(is.numeric(opt$var_filter) & opt$var_filter>0){
    is_var_filter<<-TRUE
    var_filter=opt$var_filter
  }else{
    is_var_filter<<-FALSE
    var_filter=opt$var_filter
  }
  start_genes<-ncol(gcopynumber)
    put("optional files read in",console = verbose)
    put(paste("for Species ID",s_id),console = verbose)
    put(paste("filtering with samples per copynumber",samples_per_copynumber),console = verbose)
    put(paste("depth at",depth_cutoff),console = verbose)
    put(paste("GRM used:",GRM_used),console = verbose)
    put(paste("pangenome used:",pangenome_used),console = verbose)
    put(paste("metadata used:",metadata_used),console = verbose)
    put(paste( "genes summary used:",genes_summary_used),console = verbose)
    put(paste("metadata used:",metadata_used),console = verbose)
    put(paste("log scale:",log_scale),console = verbose)
    put(paste("mean center:",mean_center),console = verbose)
    put(paste("using var filter:",is_var_filter),console = verbose)
    put(paste("starting with:", start_genes,"samples and genes:", length(unique(gcopynumber$gene_id))),console = verbose)
    
  
  prep_genes_run<-list(gcopynumber,gdepth,depth_cutoff,
                       samples_per_copynumber,verbose,opt$make_plots,opt$write_csv,
                       output_dir,s_id,pangenome_used,
                       centroid_prevalence_file,centroid_prevalence_cutoff,genes_info,GRM_used,
                       GRM,genes_summary_used,genes_summary,
                        metadata_used,metadata,min_control_case_ratio,log_scale,mean_center,is_var_filter,var_filter)
    

  
  return(prep_genes_run)
}




prep_genes_function_R<-function(gcopynumber,gdepth,depth_cutoff,samples_per_copynumber,
                              verbose=FALSE,make_plots=FALSE,write_csv=FALSE,output_dir=NULL,s_id="s_id",
                              pangenome_used=FALSE,centroid_prevalence_file=NULL,centroid_prevalence_cutoff=.7,genes_info=NULL,
                              GRM_used=FALSE,GRM=NULL,
                              genes_summary_used=FALSE,genes_summary=NULL,
                              metadata_used=FALSE,metadata=NULL,min_control_case_ratio=0,log_scale=FALSE,mean_center=FALSE,is_var_filter=FALSE,var_filter=0){
  
  list_of_samples <- colnames(gcopynumber)[-1]
  start_samples<-ncol(gcopynumber)
  start_genes<-length(unique(gcopynumber$gene_id))
  if(GRM_used){
    list_of_samples<-list_of_samples[which(list_of_samples %in% colnames(GRM))]
    if(verbose){
      put(paste("number of samples after GRM filter:",length(list_of_samples)),console = verbose)
    }
  }
  if(genes_summary_used){
    list_of_samples_gs <- genes_summary %>% filter(marker_coverage > 0) %>% filter(species_id == s_id) %>% .$sample_name %>% unique()
    list_of_samples<-list_of_samples[which(list_of_samples %in% list_of_samples_gs)]
    if(verbose){
      put(paste("number of samples after genes summary filter:",length(list_of_samples)),console = verbose)
    }
  }
  stopifnot(length(list_of_samples)>0)
  ##filter by samples first order of filtering... 1) GRM, 2) genes_summary, 3) pangeome, 4) depth and number of samples per, 5) metadata
  
  
  
  
  gdepth %<>% pivot_longer(setdiff(colnames(gdepth), "gene_id"),names_to="sample_name", values_to="gene_depth")
  gdepth %<>% filter(sample_name %in% list_of_samples)
  gcopynumber %<>% pivot_longer(setdiff(colnames(gcopynumber), "gene_id"),names_to="sample_name", values_to="copy_number")
  gcopynumber %<>% filter(sample_name %in% list_of_samples) %>% filter(copy_number>0)
  gdepth %<>% filter(gene_depth >= depth_cutoff) 
  put(paste("Gene-level first filter: gene depth >=",depth_cutoff),console = verbose)
  # compute gene occurrence frequency
  total_sample_counts<-length(list_of_samples)
  byGene <- gdepth %>% group_by(gene_id) %>% summarize(sample_counts = n()) %>% ungroup() %>%
    mutate(sample_freq = round(sample_counts / total_sample_counts, 2))
    put(paste("Filter out gene is present in less than", samples_per_copynumber, "samples (for modeling purpose)"),console = verbose)
  
  byGene %<>% filter(sample_counts >= samples_per_copynumber)
  if(pangenome_used){
    byGene<-left_join(byGene,genes_info,by=c("gene_id"="centroid_80"))
    byGene %<>% mutate(isCore = ifelse(centroid_prevalence >= centroid_prevalence_cutoff, "core", "accessory"))
    if(make_plots){
      byGene %>% ggplot(aes(x=sample_freq,y=centroid_prevalence))+  stat_bin_hex(bins=50)+labs(title=paste("For Species",s_id,"sample_freq verse centriod prevalence"))
      ggsave(file.path(output_dir, paste0(s_id,".gene_2d_histogram.pdf")), width = 7, height = 6)
      byGene %>% ggplot(aes(x = sample_freq,color=isCore)) + geom_histogram(bins = 30) + ggtitle(paste("Gene occurrence for species:", s_id))
      ggsave(file.path(output_dir, paste0(s_id,".gene_histogram.pdf")), width = 7, height = 6)
      
    }
      byGene<-byGene %>% filter(isCore=="accessory") #%>% filter(!(sample_freq>.9  & centroid_prevalence < .1))
      put(paste("labeled genes core if >=", centroid_prevalence_cutoff, "of examples for species had the gene"),console = verbose)
      put(paste("number of genes length after pangenome filter:",length(unique(byGene$gene_id))),console = verbose)
      
      }else{
    if(make_plots){
      byGene %>% ggplot(aes(x = sample_freq)) + geom_histogram(bins = 30) + ggtitle(paste("Gene occurrence for species:", s_id))
      ggsave(file.path(output_dir, paste0(s_id,".gene_histogram.pdf")), width = 7, height = 6)
    }
    
  }
  ## depth and copy number togetehr
  df <- inner_join(gdepth, gcopynumber, by=c("gene_id", "sample_name"))
  ## filter to genes with enough samples and enough depth
  df %<>% filter(gene_id %in% unique(byGene$gene_id))
  ## filter to genes with enough samples and enough depth
  df %<>% left_join(byGene, by=c("gene_id"))
  df<-df %>% mutate(copy_number=ifelse(copy_number<.4,0,copy_number))
  list_of_genes<-unique(byGene$gene_id)
  if(metadata_used){
    model_df_input<-df %>% left_join(metadata, by=c("sample_name"))
    keep_genes<-model_df_input %>% group_by(gene_id,y) %>% 
      summarize(num_samples=n(),.groups="drop") %>% 
      pivot_wider(names_from=y,values_from=num_samples, values_fill=NA)%>% mutate(control_case_ratio=`0`/`1`)
   keep_genes<-keep_genes %>% 
      filter(control_case_ratio>= min_control_case_ratio) #filter for enough samples in control and not control
    list_of_genes<-list_of_genes[which(list_of_genes %in% keep_genes$gene_id)]
      put(paste("number of genes length after metadata filter:",length(list_of_genes)),console = verbose)
  }
  copy_number_for_model<-df %>% filter(gene_id %in% list_of_genes)
  n_0<-copy_number_for_model %>% group_by(gene_id) %>% summarize(n_0=sum(copy_number==0))
  copy_number_for_model<-left_join(copy_number_for_model,n_0) %>% mutate(per_0=n_0/sample_counts) %>% filter(per_0<.3)
  if(log_scale){
    put("log_scaling",console = verbose)
    copy_number_for_model <- copy_number_for_model %>%mutate(base_copy_number=copy_number) %>%  mutate(copy_number=log(copy_number+1))
  }
  if(mean_center){
    put("mean centering",console = verbose)
    copy_number_for_model <- copy_number_for_model %>% group_by(gene_id) %>% mutate(gene_mean=mean(copy_number)) %>% mutate(copy_number=copy_number-gene_mean) 
  }
  if(is_var_filter){
    copy_number_for_model <- copy_number_for_model %>% group_by(gene_id) %>% mutate(gene_var=var(copy_number)) %>% filter(gene_var>=var_filter) 
    put(paste("number of genes left after var filter:",length(unique(copy_number_for_model$gene_id)),console = verbose))
  }
 
    put("gene filtering complete",console = verbose)
    put(paste("for Species ID",s_id),console = verbose)
    put(paste("filtering with samples per copynumber",samples_per_copynumber),console = verbose)
    put(paste("depth at",depth_cutoff),console = verbose)
    put(paste("GRM used:",GRM_used),console = verbose)
    put(paste("pangenome used:",pangenome_used),console = verbose)
    put(paste("metadata used:",metadata_used),console = verbose)
    put(paste( "genes summary used:",genes_summary_used),console = verbose)
    put(paste("metadata used:",metadata_used),console = verbose)
    put(paste("starting with:", start_genes,"genes and samples:",start_samples),console = verbose)
    put(paste("log scale:",log_scale),console = verbose)
    put(paste("mean center:",mean_center),console = verbose)
    put(paste("var filter:",is_var_filter),console = verbose)
    put(paste("var filter level:",var_filter),console = verbose)
    put(paste("ending with:", length(unique(copy_number_for_model$sample_name)),"samples"),console = verbose)
    put(paste("and with:", length(unique(copy_number_for_model$gene_id)),"genes"),console = verbose)
          
  
  df$s_id<-s_id
  if(write_csv){
    write.csv(copy_number_for_model,file.path(output_dir,paste0(s_id,".copynumber_data.csv")))
  }else{
    return(copy_number_for_model)
  }

  
}




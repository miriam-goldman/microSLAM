manhattanFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
    float dist_com=0;
    int n_sites=0;
    for (int i = 0; i < A.size(); i++){
      if(A[i]>=0 & B[i] >=0){
        dist_com += fabs(A[i]-B[i]);
        n_sites++;
      }
	          }
	          return dist_com/n_sites;
  }", depends = c("RcppArmadillo"),rebuild = TRUE)
alleleFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
    double allele_match=0;
    int n_sites=0;
    for (int i = 0; i < A.size(); i++){
      if(A[i]>=0 & B[i] >=0){
          n_sites++;
          if(fabs(A[i]-B[i])>.6){
            allele_match++;
          }
      }
	   }
	          return allele_match/n_sites;
  }", depends = c("RcppArmadillo"))




jaccardFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
    double num_a=arma::accu(A>=0);
    double num_b=arma::accu(B>=0);
    double num_both=arma::accu((B>=0 && A>=0));
    double dist;
    dist=1-(num_both/(num_a+num_b-num_both));
    return dist;
  }", depends = c("RcppArmadillo"))

#' validate_snps_input
#' 
#' helper function to validate the input from command line for snps
#' 
#' @param opt the inputs from the command line
#' @export
validate_snps_input<-function(opt){
  s_id=opt$species_id
  verbose=opt$verbose
  put(paste("species id for out files is:",s_id),console = verbose)
  if(isTRUE(!is.na(opt$snps_info_file))){
    if(file_test("-f",opt$snps_info_file)){
      snp_info<-fread(opt$snps_info_file)
    }else{
      put("snp info file doesnt exist",console = verbose)
      stop()
    }
  }else{
    put("snp info file doesnt exist",console = verbose)
    stop()
  }
  if(isTRUE(!is.na(opt$snps_freq_file))){
    if(file_test("-f",opt$snps_freq_file)){
      snp_freq<-fread(opt$snps_freq_file)
    }else{
      put("snp freq file doesnt exist",console = verbose)
     stop()
    }
  }else{
    put("snp freq file doesnt exist",console = verbose)
    stop()
  }
  
  if(isTRUE(!is.na(opt$snps_depth_file))){
    if(file_test("-f",opt$snps_depth_file)){
      snp_depth<-fread(opt$snps_depth_file)
    }else{
      put("snp depth file doesnt exist",console = verbose)
      stop()
    }
  }else{
    put("snp depth file doesnt exist",console = verbose)
    stop()
  }
  
  pangenome_used<-isTRUE(!is.na(opt$centroid_to_repgenes_file) && !is.na(opt$centroid_prevalence_file))
  
  if(pangenome_used){
      put("running pangeome step",console = verbose)
    if(file_test("-f",opt$centroid_to_repgenes_file)){
      centroid_to_repgenes<-fread(opt$centroid_to_repgenes_file)
      colnames(centroid_to_repgenes)<-c("rep_gene_id","centriod_gene_id")
    }else{
      put("centriod to rep gene doesnt exist",console = verbose)
    }
    if(file_test("-f",opt$centroid_prevalence_file)){
      core_label<-fread(opt$centroid_prevalence_file)
      colnames(core_label)<-c("centriod_gene_id","centroid_prevalence",	"centroid_ength")
    }else{
      put("centriod prevalence doesnt exist",console = verbose)
    }
    ref_core_labeled<-left_join(centroid_to_repgenes,core_label,by=c("centriod_gene_id"))
    snp_info<-snp_info %>% left_join(ref_core_labeled,by=c("gene_id"="rep_gene_id"))
    centroid_prevalence_cutoff<-opt$centroid_prevalence_cutoff
    snp_info$core<-snp_info$centroid_prevalence>centroid_prevalence_cutoff
    put("core labled with centriod prevalnce",console = verbose)
  }else{
    centroid_prevalence_cutoff<-0
    put("not running pangenome",console = verbose)
    snp_info$core<-snp_info$locus_type=="CDS"
    put("using snps in coding region",console = verbose)
  }
  
  if(isTRUE(!is.na(opt$genes_summary))){
    if(file_test("-f",opt$genes_summary)){
      if(verbose){
        put("reading in genes summary",console = verbose)
      }
      genes_summary<<-fread(opt$genes_summary) %>% filter(marker_coverage > 0)
      genes_summary_used=TRUE
    }else{
      genes_summary<-NULL
      genes_summary_used=FALSE
      put("genes depth file doesnt exist",console = verbose)
    }
  }else{
    genes_summary<-NULL
    genes_summary_used=FALSE
    put("genes depth file doesnt exist",console = verbose)
  }
  
  
  if(is.numeric(opt$abosulte_filter)){
    a=opt$abosulte_filter
  }
  if(is.numeric(opt$number_of_samples_for_sites)){
    number_of_samples_for_sites=opt$number_of_samples_for_sites
  }
  if(is.numeric(opt$sample_median_depth_filter)){
    sample_median_depth_filter=opt$sample_median_depth_filter
  }
  run_qp=opt$run_qp
  if(run_qp){
    if(is.numeric(opt$median_lower_filter)){
      l=opt$median_lower_filter
    }
    if(is.numeric(opt$median_upper_filter)){
      u=opt$median_upper_filter
    }
  }else{
    l=NA
    u=NA
  }
  
  put("optional files read in",console = verbose)
  put(paste("for Species ID",s_id),console = verbose)
  put(paste("filtering with samples per sites",number_of_samples_for_sites),console = verbose)
  put(paste("abulose depth at",a),console = verbose)
  put(paste("sample median depth filter",sample_median_depth_filter),console = verbose)
  put(paste("pangenome used:",pangenome_used),console = verbose)
  put(paste("genes summary used:",genes_summary_used),console = verbose)
  if(pangenome_used){
    put(paste("cut off for core genes",centroid_prevalence_cutoff),console = verbose)
    
  }
  if(run_qp){
    put(paste("lower median filter",l),console = verbose)
    put(paste("upper median filter",u),console = verbose)
  }
  start_snps=length(unique(snp_info$site_id))
  start_samples=ncol(snp_freq)-1
  put(paste("starting with:", start_samples,"samples and snps:",start_snps),console = verbose)
  
  snp_info %>% dplyr::count(snp_type, locus_type, site_type) %>% 
    spread(site_type, n, fill = 0) %>% put(console = verbose)
  
  
  make_plots=opt$make_plots
  prep_snps_inputs<-list(snp_freq,snp_depth,snp_info,sample_median_depth_filter,number_of_samples_for_sites,verbose,make_plots,pangenome_used,centroid_prevalence_cutoff,run_qp,l,u,a,start_samples,start_snps,genes_summary_used,genes_summary)
  return(prep_snps_inputs)
}

#' prep_snps_function_R
#' 
#' helper function to filter snps data from midas2 and make GRM
#' 
#' @param snp_freq snp_freq from MIDAS2
#' @param snp_depth snp_depth from MIDAS2
#' @param snp_info snp_info from MIDAS2
#' @param sample_median_depth_filter sample median depth filter
#' @param number_of_samples_for_sites number of samples for sites
#' @export
prep_snps_function_R<-function(snp_freq,snp_depth,snp_info,sample_median_depth_filter,number_of_samples_for_sites,verbose=FALSE,make_plots=FALSE,pangenome_used=FALSE,centroid_prevalence_cutoff=.8,run_qp=FALSE,l=.3,u=3,a=5,start_samples,start_snps,genes_summary_used,genes_summary){
  ######## Filter Samples based on D_median_cds 
  # Compute median site depth for all protein coding genes
  manhattanFuncPtr <- cppXPtr(
    "double customDist(const arma::mat &A, const arma::mat &B) {
    float dist_com=0;
    int n_sites=0;
    for (int i = 0; i < A.size(); i++){
      if(A[i]>=0 & B[i] >=0){
        dist_com += fabs(A[i]-B[i]);
        n_sites++;
      }
	          }
	          return dist_com/n_sites;
  }", depends = c("RcppArmadillo"),rebuild = TRUE)
  snp_depth <- snp_depth %>%
    filter(site_id %in% unique(snp_info$site_id))
  nsamples = ncol(snp_depth)-1
  list_of_samples <- colnames(snp_depth)[-1]
  put(paste("number of samples",nsamples),console = verbose)
  # Site not covered in a certain <sample, species> pair is site-depth = 0, and it should not be included in the computation.
  snp_depth[snp_depth == 0] <- NA
  per_sample_median_depth <- apply(snp_depth[,-1], 2, function(x) median(x, na.rm =T))
  samples_pass_depth <- names(per_sample_median_depth[per_sample_median_depth >= sample_median_depth_filter]) #<-----
  list_of_samples<-list_of_samples[which(list_of_samples %in% samples_pass_depth)]
  snp_depth %<>% select(site_id, matches(list_of_samples))
  if(genes_summary_used){
    list_of_samples_gs <- genes_summary %>% filter(species_id == s_id) %>% .$sample_name %>% unique()
    list_of_samples<-list_of_samples[which(list_of_samples %in% list_of_samples_gs)]
    snp_depth %<>% select(site_id, matches(list_of_samples))
  }
  nsamples2 = ncol(snp_depth)-1
  put(paste("number of samples after filter",nsamples2),console = verbose)
  D <- data.frame(sample_name = names(per_sample_median_depth), median_site_depth=per_sample_median_depth, row.names=NULL)
  if(make_plots){
    g2 <- D %>% 
      ggplot(aes(x = median_site_depth)) + geom_histogram(bins=30) + geom_vline(xintercept = sample_median_depth_filter, color = "red") + 
      scale_x_log10() + 
      ggtitle(paste(s_id, ":", nsamples2, "out of", nsamples, "samples passing the filter at:",sample_median_depth_filter))
      put("number of samples passing filter plot",console=verbose)
    
    ggsave(file.path(output_dir, paste0(s_id,".depth_histogram.pdf")),g2, width = 7, height = 6)
    
  }
 
  
  
  ######### Keep sites [0.3, 3] * D_median_cds
  
  
  
  if(run_qp){
    info_for_qp <- snp_info %>% filter(site_type == "4D")
    depth_for_qp<-snp_depth %>% filter(site_id %in% unique(info_for_qp$site_id)) 
    depth_for_qp %<>% pivot_longer(matches(list_of_samples),names_to="sample_name", values_to="site_depth")
    depth_for_qp %<>% 
      left_join(D %>% select(sample_name, median_site_depth),by=c("sample_name")) %>%
      mutate(min_bound =  l* median_site_depth, max_bound = u * median_site_depth)
    depth_for_qp %<>%
      filter(site_depth >= min_bound & site_depth <= max_bound) %>%
      filter(site_depth >= a) %>% #<-------------
    select(site_id, sample_name, site_depth, median_site_depth)
    sc_df <- depth_for_qp %>%
      group_by(site_id) %>% 
      summarise(sample_counts = n(),total_reads=sum(site_depth)) %>%
      ungroup()
    sc_df$pos<-str_split_i(sc_df$site_id,"\\|",4)
    sc_df$contig<-str_split_i(sc_df$site_id,"\\|",3)
    sites_no <- unique(sc_df %>% filter(sample_counts <= number_of_samples_for_sites) %>% .$site_id)
    stopifnot(nrow(sites_no) == 0)
    if (length(sites_no) > 0) {
      depth_for_qp %<>% filter(!site_id %in% sites_no)
    }
    freq_for_qp <- snp_freq%>%
      filter(site_id %in% unique(depth_for_qp$site_id)) %>%
      select(site_id, matches(samples_pass_depth))
    ######### Read in population minor allele frequency
    
    
    freq_for_qp %<>%
      pivot_longer(matches(samples_pass_depth),names_to="sample_name", values_to="allele_freq") %>%
      filter(allele_freq != -1)
    df <- left_join(depth_for_qp, freq_for_qp, by=c("site_id", "sample_name"))
    df %<>% 
      mutate(allele_direction = ifelse(allele_freq <= 0.2, "low", "int")) %>%
      mutate(allele_direction = ifelse(allele_freq >= 0.8, "high", allele_direction))
    
    
    f2 <- file.path(output_dir, paste0(s_id,".site_depth_with_allele_freq.pdf"))
    t1 <- df %>% select(site_depth, allele_freq) %>%
      group_by(site_depth, allele_freq) %>%
      dplyr::count() %>%
      dplyr::rename(occurence = n) %>%
      ungroup()
    ######### Compute ND and QP
    by_site <- df %>% 
      group_by(site_id) %>%
      summarise(n_low = sum(allele_direction == "low"), n_high = sum(allele_direction == "high")) %>%
      ungroup() %>%
      mutate(voted_direction = ifelse(n_low > n_high, "low", "high")) %>%
      mutate(voted_counted = ifelse(n_low > n_high, n_low, n_high))
    by_site %<>% left_join(sc_df, by=c("site_id")) %>%
      mutate(voted_freq = round(voted_counted / sample_counts, 4)) %>%
      select(-voted_counted) %>%
      select(site_id, sample_counts, everything())
    gene_info_1<-inner_join(by_site,info_for_qp,by=c("site_id"))
      put(print(head(info_for_qp)),console = verbose)
      put(print(head(gene_info_1)),console = verbose)
    if(make_plots & pangenome_used){
      gene_info_1<-gene_info_1 %>% group_by(gene_id) %>% summarize(mean_sample_count=mean(sample_counts.x),sum_total_reads=sum(total_reads),mean_total_reads=mean(total_reads), centroid_prevalence=mean(centroid_prevalence))
      
      gene_plot<-gene_info_1 %>% ggplot(aes(mean_sample_count,centroid_prevalence,color=mean_total_reads))+geom_point()+labs(title=paste("mean sample count for gene verse centroid prevalence for species",s_id))
      gene_plot_1<-ggExtra::ggMarginal(gene_plot, type = "histogram")
      print(gene_plot_1,newpage = TRUE)
      f5 <- file.path(output_dir, paste(s_id,"gene_hist.pdf", sep=""))
      ggsave(f5, gene_plot_1, width = 7, height = 5)
    }
   
    by_sample <- df %>%
      group_by(sample_name) %>%
      summarise(n_int = sum(allele_direction == "int"), n_total = n()) %>%
      ungroup()
    
    denominator  <- df %>%
      filter(allele_direction != "int") %>%
      left_join(by_site %>% select(site_id, voted_direction, voted_freq), by=c("site_id")) %>%
      mutate(panel_freq = ifelse(allele_direction == voted_direction, voted_freq, 1-voted_freq))
    
    nd <- denominator %>% 
      group_by(sample_name) %>%
      summarise(Nd = sum(panel_freq)) %>%
      ungroup()
    
    by_sample %<>% left_join(nd, by=c("sample_name")) %>%
      mutate(QP = round(n_int / Nd, 3))
    if(make_plots){
      f3 <- file.path(output_dir, paste(s_id,"QP_hist.pdf", sep=""))
      g3 <- by_sample %>%
        ggplot(aes(x = QP)) + geom_histogram(bins=30) +
        theme_bw() + geom_vline(xintercept = 0.1, color = "red")
      print(g3)
      ggsave(f3, g3, width = 7, height = 5)
    }
    
    put(paste("printing to",file.path(output_dir, paste0(s_id,"by_sample.tsv"))),console = verbose)
    by_sample %>%
      write.table(file.path(output_dir, paste0(s_id,"by_sample.tsv")), sep = "\t", quote = F, row.names = F)
      put("QP completed",console = verbose)
      put(head(by_sample),console = verbose)
  }
  
  #### re-filter everything so we arent just looking at 4D sites
  
  info_for_distance<-snp_info %>% filter(core==TRUE) %>% filter(snp_type=="bi")
  depth_for_distance<-snp_depth %>% filter(site_id %in% unique(info_for_distance$site_id)) 
  depth_for_distance %<>% pivot_longer(matches(list_of_samples), names_to="sample_name", values_to="site_depth")
  
  depth_for_distance %<>% 
    left_join(D %>% select(sample_name, median_site_depth),by=c("sample_name"))
  
  depth_for_distance %<>%
    filter(site_depth >= a) %>% #<-------------
  select(site_id, sample_name, site_depth, median_site_depth)
  
  
  sc_df_2 <- depth_for_distance %>%
    group_by(site_id) %>% 
    summarise(sample_counts = n(),total_reads=sum(site_depth)) %>%
    ungroup()
  sc_df_2$pos<-str_split_i(sc_df_2$site_id,"\\|",4)
  sc_df_2$contig<-str_split_i(sc_df_2$site_id,"\\|",3)
  if(make_plots){
    put("sample depth historgram",console = verbose)
    sample_counts_per_site<-sc_df_2 %>% ggplot(aes(x = sample_counts)) + geom_histogram(bins=30)+labs(title="sample depth histogram for sites")
    
    f4 <- file.path(output_dir, paste(s_id,".site_hist.pdf", sep=""))
    ggsave(f4, sample_counts_per_site, width = 7, height = 5)
    
  }
 
  sites_pass <- unique(sc_df_2 %>% filter(sample_counts >= number_of_samples_for_sites) %>% .$site_id)
  stopifnot(nrow(sites_pass) > 0)
  if (length(sites_pass) > 0) {
    depth_for_distance %<>% filter(site_id %in% sites_pass)
  }
  ######### Read in population minor allele frequency
  freq_for_distance <- snp_freq%>%
    filter(site_id %in% unique(depth_for_distance$site_id)) %>%
    select(site_id, matches(samples_pass_depth))
  
  freq_for_distance %<>%
    pivot_longer(matches(samples_pass_depth),names_to="sample_name", values_to="allele_freq") %>% filter(allele_freq > -1)
  
  
  df_for_distance <- left_join(depth_for_distance, freq_for_distance, by=c("site_id", "sample_name"))
  site_df<-df_for_distance %>% select(site_id,sample_name,allele_freq) %>% pivot_wider(names_from = site_id,values_from=allele_freq)
  put(head(site_df),console = verbose)
  df_for_distance  %>%  write.table(file.path(output_dir, paste0(s_id,".snp_sample_site.tsv")), sep = "\t", quote = F,row.names=FALSE)
  freq_mat_dist_man<-parDist(as.matrix(site_df[,-1]), method="custom", func = manhattanFuncPtr)
  freq_mat_dist_man<-as.matrix(freq_mat_dist_man)
  
  dimnames(freq_mat_dist_man)<-c(site_df[,1],site_df[,1])
  freq_mat_GRM_man<-1-freq_mat_dist_man
  freq_mat_GRM_man  %>%  write.table(file.path(output_dir, paste0(s_id,".GRM.tsv")), sep = "\t", quote = F)
  freq_mat_GRM_man  %>%  write.table(file.path(output_dir, paste0(s_id,".distance.tsv")), sep = "\t", quote = F)
  if(make_plots){
    heatmap_file_name<-file.path(output_dir, paste0(s_id,".heatmap.pdf"))
    if(run_qp){
      annotation_col = data.frame(qp=as.character(by_sample$QP<.1))
      rownames(annotation_col) = by_sample$sample_name
      ann_colors = list(
        qp = c("TRUE" = "#1B9E77", "FALSE" = "#D95F02")
      )
      
      pheatmap(freq_mat_dist_man,annotation_row = annotation_col,annotation_col=annotation_col,annotation_colors = ann_colors,filename=heatmap_file_name)
      
    }else{
      pheatmap(freq_mat_dist_man,filename=heatmap_file_name)
      
    }
  }
  end_snps=ncol(site_df)-1
  end_samples=nrow(site_df)
  
  
  put(paste("for Species ID",s_id),console = verbose)
  put(paste("filtering with samples per sites",number_of_samples_for_sites),console = verbose)
  put(paste("abulose depth at",a),console = verbose)
  put(paste("sample median depth filter",sample_median_depth_filter),console = verbose)
  put(paste("pangenome used:",pangenome_used),console = verbose)
  if(pangenome_used){
    put(paste("cut off for core genes",centroid_prevalence_cutoff),console = verbose)
    
  }
  if(run_qp){
    put(paste("lower median filter",l),console = verbose)
    put(paste("upper median filter",u),console = verbose)
  }
  start_snps=length(unique(snp_info$site_id))
  start_samples=ncol(snp_freq)-1
  put(paste("starting with:", start_samples,"samples and snps:",start_snps),console = verbose)
  put(paste("ending with:", end_samples,"samples and snps:",end_snps),console = verbose)
  
 
}

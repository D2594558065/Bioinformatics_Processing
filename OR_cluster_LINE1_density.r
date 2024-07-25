calculate_line1_density <- function(human_path, line1_path, distance_threshold = 1000000) {
  library(dplyr)

  # 读取和处理人类数据
  dfh <- read.csv(human_path, sep = "\t") %>% 
    filter(grepl("olfactory receptor", description)) %>% 
    dplyr::select(Symbol, genomic_nucleotide_accession.version, start_position_on_the_genomic_accession, end_position_on_the_genomic_accession)
  
  colnames(dfh) <- c("Symbol", "chromosome", "start_position_on_the_genomic_accession", "end_position_on_the_genomic_accession")
  dfh <- na.omit(dfh)
  dfh <- dfh[order(dfh$chromosome, dfh$start_position_on_the_genomic_accession),]

  # 初始化簇ID
  cluster_id <- 1
  dfh$cluster <- NA
  dfh$cluster[1] <- cluster_id

  # 聚类逻辑
  for (i in 2:nrow(dfh)) {
    if (dfh$chromosome[i] == dfh$chromosome[i - 1] && 
        (dfh$start_position_on_the_genomic_accession[i] - dfh$end_position_on_the_genomic_accession[i - 1]) <= distance_threshold) {
      dfh$cluster[i] <- cluster_id
    } else {
      cluster_id <- cluster_id + 1
      dfh$cluster[i] <- cluster_id
    }
  }

  dfh <- dfh %>%
    group_by(cluster) %>%           
    filter(n() >= 5) %>%            
    ungroup()

  cluster_boundaries <- dfh %>%
    group_by(cluster, chromosome) %>%  
    summarise(
      min_start = min(start_position_on_the_genomic_accession),
      max_end = max(end_position_on_the_genomic_accession),
      .groups = 'drop'
    )

  # 读取和处理LINE1数据
  line1_bed <- read.csv(line1_path, sep = "", skip = 2, header = FALSE) 

  line1_bed <- line1_bed %>% 
  	filter(grepl("LINE", V11)) %>% 
    	dplyr::select(V5, V6, V7, V10, V11)


  colnames(line1_bed) <- c("chromosome", "start", "end", "family", "superfamily")

  # 计算每个簇中的LINE1总长度和密度
  cluster_line1_density <- cluster_boundaries %>%
    rowwise() %>%
    mutate(
      line1_length = sum(
        pmax(0, pmin(line1_bed$end[line1_bed$chromosome == chromosome & 
                                      line1_bed$start <= max_end & 
                                      line1_bed$end >= min_start], max_end) - 
              pmax(line1_bed$start[line1_bed$chromosome == chromosome & 
                                           line1_bed$start <= max_end & 
                                           line1_bed$end >= min_start], min_start) + 1)
      ),
      cluster_length = max_end - min_start + 1,  
      density = line1_length / cluster_length  
    ) %>%
    ungroup()

  return (cluster_line1_density)
}


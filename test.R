bat <- function(analysis = "paired", tumourname, normalname, tumour_data_file, 
                normal_data_file, imputeinfofile, g1000prefix, problemloci, 
                gccorrectprefix = NULL, repliccorrectprefix = NULL, g1000allelesprefix = NA, 
                ismale = NA, data_type = "wgs", impute_exe = "impute2", 
                allelecounter_exe = "alleleCounter", nthreads = 8, platform_gamma = 1, 
                phasing_gamma = 1, segmentation_gamma = 10, segmentation_kmin = 3, 
                phasing_kmin = 1, clonality_dist_metric = 0, ascat_dist_metric = 1, 
                min_ploidy = 1.6, max_ploidy = 4.8, min_rho = 0.1, min_goodness = 0.63, 
                uninformative_BAF_threshold = 0.51, min_normal_depth = 10, 
                min_base_qual = 20, min_map_qual = 35, calc_seg_baf_option = 3, 
                skip_allele_counting = F, skip_preprocessing = F, skip_phasing = F, 
                externalhaplotypefile = NA, usebeagle = FALSE, beaglejar = NA, 
                beagleref.template = NA, beagleplink.template = NA, beaglemaxmem = 10, 
                beaglenthreads = 1, beaglewindow = 40, beagleoverlap = 4, 
                javajre = "java", write_battenberg_phasing = T, multisample_relative_weight_balanced = 0.25, 
                multisample_maxlag = 100, segmentation_gamma_multisample = 5, 
                snp6_reference_info_file = NA, apt.probeset.genotype.exe = "apt-probeset-genotype", 
                apt.probeset.summarize.exe = "apt-probeset-summarize", norm.geno.clust.exe = "normalize_affy_geno_cluster.pl", 
                birdseed_report_file = "birdseed.report.txt", heterozygousFilter = "none", 
                prior_breakpoints_file = NULL, GENOMEBUILD = "hg19", chrom_coord_file = NULL){
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  if (analysis == "cell_line") {
    calc_seg_baf_option = 1
    phasing_gamma = 1
    phasing_kmin = 2
    segmentation_gamma = 20
    segmentation_kmin = 3
    normalname = paste0(tumourname, "_normal")
  }
  if (analysis == "germline") {
    calc_seg_baf_option = 1
    phasing_gamma = 3
    phasing_kmin = 1
    segmentation_gamma = 3
    segmentation_kmin = 3
  }
  if (data_type == "wgs" & is.na(ismale)) {
    stop("Please provide a boolean denominator whether this sample represents a male donor")
  }
  if (data_type == "wgs" & is.na(g1000allelesprefix)) {
    stop("Please provide a path to 1000 Genomes allele reference files")
  }
  if (data_type == "wgs" & is.null(gccorrectprefix)) {
    stop("Please provide a path to GC content reference files")
  }
  if (!file.exists(problemloci)) {
    stop("Please provide a path to a problematic loci file")
  }
  if (!file.exists(imputeinfofile)) {
    stop("Please provide a path to an impute info file")
  }
  check.imputeinfofile(imputeinfofile = imputeinfofile, is.male = ismale, 
                       usebeagle = usebeagle)
  nsamples <- length(tumourname)
  if (nsamples > 1) {
    if (length(skip_allele_counting) < nsamples) {
      skip_allele_counting = rep(skip_allele_counting[1], 
                                 nsamples)
    }
    if (length(skip_preprocessing) < nsamples) {
      skip_preprocessing = rep(skip_preprocessing[1], 
                               nsamples)
    }
    if (length(skip_phasing) < nsamples) {
      skip_phasing = rep(skip_phasing[1], nsamples)
    }
  }
  if (data_type == "wgs" | data_type == "WGS") {
    if (nsamples > 1) {
      print(paste0("Running Battenberg in multisample mode on ", 
                   nsamples, " samples: ", paste0(tumourname, collapse = ", ")))
    }
    chrom_names = get.chrom.names(imputeinfofile, ismale, 
                                  analysis = analysis)
  }
  else if (data_type == "snp6" | data_type == "SNP6") {
    if (nsamples > 1) {
      stop(paste0("Battenberg multisample mode has not been tested with SNP6 data"))
    }
    chrom_names = get.chrom.names(imputeinfofile, TRUE)
    logr_file = paste(tumourname, "_mutantLogR.tab", sep = "")
    allelecounts_file = NULL
  }
  print(chrom_names)
  for (sampleidx in 1:nsamples) {
    if (!skip_preprocessing[sampleidx]) {
      if (data_type == "wgs" | data_type == "WGS") {
        clp = parallel::makeCluster(nthreads)
        doParallel::registerDoParallel(clp)
        if (analysis == "paired") {
          prepare_wgs(chrom_names = chrom_names, tumourbam = tumour_data_file[sampleidx], 
                      normalbam = normal_data_file, tumourname = tumourname[sampleidx], 
                      normalname = normalname, g1000allelesprefix = g1000allelesprefix, 
                      g1000prefix = g1000prefix, gccorrectprefix = gccorrectprefix, 
                      repliccorrectprefix = repliccorrectprefix, 
                      min_base_qual = min_base_qual, min_map_qual = min_map_qual, 
                      allelecounter_exe = allelecounter_exe, min_normal_depth = min_normal_depth, 
                      nthreads = nthreads, skip_allele_counting = skip_allele_counting[sampleidx], 
                      skip_allele_counting_normal = (sampleidx > 
                                                       1))
        }
        else if (analysis == "cell_line") {
          prepare_wgs_cell_line(chrom_names = chrom_names, 
                                chrom_coord = chrom_coord_file, tumourbam = tumour_data_file, 
                                tumourname = tumourname, g1000lociprefix = g1000prefix, 
                                g1000allelesprefix = g1000allelesprefix, 
                                gamma_ivd = 100000, kmin_ivd = 50, centromere_noise_seg_size = 1000000, 
                                centromere_dist = 500000, min_het_dist = 100000, 
                                gamma_logr = 100, length_adjacent = 50000, 
                                gccorrectprefix = gccorrectprefix, repliccorrectprefix = repliccorrectprefix, 
                                min_base_qual = min_base_qual, min_map_qual = min_map_qual, 
                                allelecounter_exe = allelecounter_exe, min_normal_depth = min_normal_depth, 
                                skip_allele_counting = skip_allele_counting[sampleidx])
        }
        else if (analysis == "germline") {
        }
        parallel::stopCluster(clp)
      }
      else if (data_type == "snp6" | data_type == "SNP6") {
        prepare_snp6(tumour_cel_file = tumour_data_file[sampleidx], 
                     normal_cel_file = normal_data_file, tumourname = tumourname[sampleidx], 
                     chrom_names = chrom_names, snp6_reference_info_file = snp6_reference_info_file, 
                     apt.probeset.genotype.exe = apt.probeset.genotype.exe, 
                     apt.probeset.summarize.exe = apt.probeset.summarize.exe, 
                     norm.geno.clust.exe = norm.geno.clust.exe, 
                     birdseed_report_file = birdseed_report_file)
      }
      else {
        print("Unknown data type provided, please provide wgs or snp6")
        q(save = "no", status = 1)
      }
    }
    if (data_type == "snp6" | data_type == "SNP6") {
      gender = infer_gender_birdseed(birdseed_report_file)
      ismale = gender == "male"
    }
}
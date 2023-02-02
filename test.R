bat <- function (analysis = "paired", tumourname, normalname, tumour_data_file, 
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
                 prior_breakpoints_file = NULL, GENOMEBUILD = "hg19", chrom_coord_file = NULL) 
{
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
  Battenberg:::check.imputeinfofile(imputeinfofile = imputeinfofile, is.male = ismale, 
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
    chrom_names = Battenberg:::get.chrom.names(imputeinfofile, ismale, 
                                  analysis = analysis)
  }
  else if (data_type == "snp6" | data_type == "SNP6") {
    if (nsamples > 1) {
      stop(paste0("Battenberg multisample mode has not been tested with SNP6 data"))
    }
    chrom_names = Battenberg:::get.chrom.names(imputeinfofile, TRUE)
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
          Battenberg:::prepare_wgs(chrom_names = chrom_names, tumourbam = tumour_data_file[sampleidx], 
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
          Battenberg:::prepare_wgs_cell_line(chrom_names = chrom_names, 
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
        Battenberg:::prepare_snp6(tumour_cel_file = tumour_data_file[sampleidx], 
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
      gender = Battenberg:::infer_gender_birdseed(birdseed_report_file)
      ismale = gender == "male"
    }
    if (!skip_phasing[sampleidx]) {
      if (!is.na(externalhaplotypefile) && file.exists(externalhaplotypefile)) {
        externalhaplotypeprefix <- paste0(normalname, 
                                          "_external_haplotypes_chr")
        if (any(!file.exists(paste0(externalhaplotypeprefix, 
                                    1:length(chrom_names), ".vcf")))) {
          print(paste0("Splitting external phasing data from ", 
                       externalhaplotypefile))
          Battenberg:::split_input_haplotypes(chrom_names = chrom_names, 
                                 externalhaplotypefile = externalhaplotypefile, 
                                 outprefix = externalhaplotypeprefix)
        }
        else {
          print("No need to split, external haplotype files per chromosome found")
        }
      }
      else {
        externalhaplotypeprefix <- NA
      }
      clp = parallel::makeCluster(nthreads)
      doParallel::registerDoParallel(clp)
      foreach::foreach(i = 1:length(chrom_names)) %dopar% 
        {
          chrom = chrom_names[i]
          print(chrom)
          Battenberg:::run_haplotyping(chrom = chrom, tumourname = tumourname[sampleidx], 
                          normalname = normalname, ismale = ismale, 
                          imputeinfofile = imputeinfofile, problemloci = problemloci, 
                          impute_exe = impute_exe, min_normal_depth = min_normal_depth, 
                          chrom_names = chrom_names, snp6_reference_info_file = snp6_reference_info_file, 
                          heterozygousFilter = heterozygousFilter, 
                          usebeagle = usebeagle, beaglejar = beaglejar, 
                          beagleref = gsub("CHROMNAME", chrom, beagleref.template), 
                          beagleplink = gsub("CHROMNAME", chrom, beagleplink.template), 
                          beaglemaxmem = beaglemaxmem, beaglenthreads = beaglenthreads, 
                          beaglewindow = beaglewindow, beagleoverlap = beagleoverlap, 
                          externalhaplotypeprefix = externalhaplotypeprefix, 
                          use_previous_imputation = (sampleidx > 1))
        }
      parallel::stopCluster(clp)
      Battenberg:::combine.baf.files(inputfile.prefix = paste(tumourname[sampleidx], 
                                                 "_chr", sep = ""), inputfile.postfix = "_heterozygousMutBAFs_haplotyped.txt", 
                        outputfile = paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", 
                                           sep = ""), chr_names = chrom_names)
    }
    Battenberg:::segment.baf.phased(samplename = tumourname[sampleidx], 
                       inputfile = paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", 
                                         sep = ""), outputfile = paste(tumourname[sampleidx], 
                                                                       ".BAFsegmented.txt", sep = ""), prior_breakpoints_file = prior_breakpoints_file, 
                       gamma = segmentation_gamma, phasegamma = phasing_gamma, 
                       kmin = segmentation_kmin, phasekmin = phasing_kmin, 
                       calc_seg_baf_option = calc_seg_baf_option)
    if (nsamples > 1 | write_battenberg_phasing) {
      Battenberg:::write_battenberg_phasing(tumourname = tumourname[sampleidx], 
                               SNPfiles = paste0(tumourname[sampleidx], "_alleleFrequencies_chr", 
                                                 chrom_names, ".txt"), imputedHaplotypeFiles = paste0(tumourname[sampleidx], 
                                                                                                      "_impute_output_chr", chrom_names, "_allHaplotypeInfo.txt"), 
                               bafsegmented_file = paste0(tumourname[sampleidx], 
                                                          ".BAFsegmented.txt"), outprefix = paste0(tumourname[sampleidx], 
                                                                                                   "_Battenberg_phased_chr"), chrom_names = chrom_names, 
                               include_homozygous = F)
    }
  }
  # if (nsamples > 1) {
  #   print("Constructing multisample phasing")
  #   multisamplehaplotypeprefix <- paste0(normalname, "_multisample_haplotypes_chr")
  #   clp = parallel::makeCluster(nthreads)
  #   doParallel::registerDoParallel(clp)
  #   foreach::foreach(i = 1:length(chrom_names)) %dopar% 
  #     {
  #       chrom = chrom_names[i]
  #       print(chrom)
  #       get_multisample_phasing(chrom = chrom, bbphasingprefixes = paste0(tumourname, 
  #                                                                         "_Battenberg_phased_chr"), maxlag = multisample_maxlag, 
  #                               relative_weight_balanced = multisample_relative_weight_balanced, 
  #                               outprefix = multisamplehaplotypeprefix)
  #     }
  #   for (sampleidx in 1:nsamples) {
  #     MutBAFfiles <- paste0(tumourname[sampleidx], "_chr", 
  #                           chrom_names, "_heterozygousMutBAFs_haplotyped.txt")
  #     heterozygousdatafiles <- paste0(tumourname[sampleidx], 
  #                                     "_chr", chrom_names, "_heterozygousData.png")
  #     raffiles <- paste0(tumourname[sampleidx], "_RAFseg_chr", 
  #                        chrom_names, ".png")
  #     segfiles <- paste0(tumourname[sampleidx], "_segment_chr", 
  #                        chrom_names, ".png")
  #     haplotypedandbafsegmentedfiles <- paste0(tumourname[sampleidx], 
  #                                              c("_heterozygousMutBAFs_haplotyped.txt", ".BAFsegmented.txt"))
  #     file.copy(from = MutBAFfiles, to = gsub(pattern = ".txt$", 
  #                                             replacement = "_noMulti.txt", x = MutBAFfiles), 
  #               overwrite = T)
  #     file.copy(from = heterozygousdatafiles, to = gsub(pattern = ".png$", 
  #                                                       replacement = "_noMulti.png", x = heterozygousdatafiles), 
  #               overwrite = T)
  #     file.copy(from = raffiles, to = gsub(pattern = ".png$", 
  #                                          replacement = "_noMulti.png", x = raffiles), 
  #               overwrite = T)
  #     file.copy(from = segfiles, to = gsub(pattern = ".png$", 
  #                                          replacement = "_noMulti.png", x = segfiles), 
  #               overwrite = T)
  #     file.copy(from = haplotypedandbafsegmentedfiles, 
  #               to = gsub(pattern = ".txt$", replacement = "_noMulti.txt", 
  #                         x = haplotypedandbafsegmentedfiles), overwrite = T)
  #     foreach::foreach(i = 1:length(chrom_names)) %dopar% 
  #       {
  #         chrom = chrom_names[i]
  #         print(chrom)
  #         input_known_haplotypes(chrom = chrom, chrom_names = chrom_names, 
  #                                imputedHaplotypeFile = paste0(tumourname[sampleidx], 
  #                                                              "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"), 
  #                                externalHaplotypeFile = paste0(multisamplehaplotypeprefix, 
  #                                                               chrom, ".vcf"), oldfilesuffix = "_noMulti.txt")
  #         GetChromosomeBAFs(chrom = chrom, SNP_file = paste0(tumourname[sampleidx], 
  #                                                            "_alleleFrequencies_chr", chrom, ".txt"), 
  #                           haplotypeFile = paste0(tumourname[sampleidx], 
  #                                                  "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"), 
  #                           samplename = tumourname[sampleidx], outfile = paste0(tumourname[sampleidx], 
  #                                                                                "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"), 
  #                           chr_names = chrom_names, minCounts = min_normal_depth)
  #         plot.haplotype.data(haplotyped.baf.file = paste0(tumourname[sampleidx], 
  #                                                          "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"), 
  #                             imageFileName = paste0(tumourname[sampleidx], 
  #                                                    "_chr", chrom, "_heterozygousData.png"), 
  #                             samplename = tumourname[sampleidx], chrom = chrom, 
  #                             chr_names = chrom_names)
  #       }
  #   }
  #   parallel::stopCluster(clp)
  #   for (sampleidx in 1:nsamples) {
  #     combine.baf.files(inputfile.prefix = paste0(tumourname[sampleidx], 
  #                                                 "_chr"), inputfile.postfix = "_heterozygousMutBAFs_haplotyped.txt", 
  #                       outputfile = paste0(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt"), 
  #                       chr_names = chrom_names)
  #   }
  #   segment.baf.phased.multisample(samplename = tumourname, 
  #                                  inputfile = paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", 
  #                                                    sep = ""), outputfile = paste(tumourname, ".BAFsegmented.txt", 
  #                                                                                  sep = ""), prior_breakpoints_file = prior_breakpoints_file, 
  #                                  gamma = segmentation_gamma_multisample, calc_seg_baf_option = calc_seg_baf_option, 
  #                                  GENOMEBUILD = GENOMEBUILD)
  # }
  # clp = parallel::makeCluster(min(nthreads, nsamples))
  # doParallel::registerDoParallel(clp)
  # foreach::foreach(sampleidx = 1:nsamples) %dopar% {
  #   print(paste0("Fitting final copy number and calling subclones for sample ", 
  #                tumourname[sampleidx]))
  #   if (data_type == "wgs" | data_type == "WGS") {
  #     logr_file = paste(tumourname[sampleidx], "_mutantLogR_gcCorrected.tab", 
  #                       sep = "")
  #     if (analysis == "paired") {
  #       allelecounts_file = paste(tumourname[sampleidx], 
  #                                 "_alleleCounts.tab", sep = "")
  #     }
  #     else {
  #       allelecounts_file = NULL
  #     }
  #   }
  #   fit.copy.number(samplename = tumourname[sampleidx], 
  #                   outputfile.prefix = paste(tumourname[sampleidx], 
  #                                             "_", sep = ""), inputfile.baf.segmented = paste(tumourname[sampleidx], 
  #                                                                                             ".BAFsegmented.txt", sep = ""), inputfile.baf = paste(tumourname[sampleidx], 
  #                                                                                                                                                   "_mutantBAF.tab", sep = ""), inputfile.logr = logr_file, 
  #                   dist_choice = clonality_dist_metric, ascat_dist_choice = ascat_dist_metric, 
  #                   min.ploidy = min_ploidy, max.ploidy = max_ploidy, 
  #                   min.rho = min_rho, min.goodness = min_goodness, 
  #                   uninformative_BAF_threshold = uninformative_BAF_threshold, 
  #                   gamma_param = platform_gamma, use_preset_rho_psi = F, 
  #                   preset_rho = NA, preset_psi = NA, read_depth = 30, 
  #                   analysis = analysis)
  #   callSubclones(sample.name = tumourname[sampleidx], baf.segmented.file = paste(tumourname[sampleidx], 
  #                                                                                 ".BAFsegmented.txt", sep = ""), logr.file = logr_file, 
  #                 rho.psi.file = paste(tumourname[sampleidx], "_rho_and_psi.txt", 
  #                                      sep = ""), output.file = paste(tumourname[sampleidx], 
  #                                                                     "_subclones.txt", sep = ""), output.figures.prefix = paste(tumourname[sampleidx], 
  #                                                                                                                                "_subclones_chr", sep = ""), output.gw.figures.prefix = paste(tumourname[sampleidx], 
  #                                                                                                                                                                                              "_BattenbergProfile", sep = ""), masking_output_file = paste(tumourname[sampleidx], 
  #                                                                                                                                                                                                                                                           "_segment_masking_details.txt", sep = ""), prior_breakpoints_file = prior_breakpoints_file, 
  #                 chr_names = chrom_names, gamma = platform_gamma, 
  #                 segmentation.gamma = NA, siglevel = 0.05, maxdist = 0.01, 
  #                 noperms = 1000, calc_seg_baf_option = calc_seg_baf_option)
  #   if (ismale & "X" %in% chrom_names) {
  #     callChrXsubclones(TUMOURNAME = tumourname[sampleidx], 
  #                       X_GAMMA = 1000, X_KMIN = 100, GENOMEBUILD = GENOMEBUILD, 
  #                       AR = TRUE)
  #   }
  #   make_posthoc_plots(samplename = tumourname[sampleidx], 
  #                      logr_file = logr_file, subclones_file = paste(tumourname[sampleidx], 
  #                                                                    "_subclones.txt", sep = ""), rho_psi_file = paste(tumourname[sampleidx], 
  #                                                                                                                      "_rho_and_psi.txt", sep = ""), bafsegmented_file = paste(tumourname[sampleidx], 
  #                                                                                                                                                                               ".BAFsegmented.txt", sep = ""), logrsegmented_file = paste(tumourname[sampleidx], 
  #                                                                                                                                                                                                                                          ".logRsegmented.txt", sep = ""), allelecounts_file = allelecounts_file)
  #   cnfit_to_refit_suggestions(samplename = tumourname[sampleidx], 
  #                              subclones_file = paste(tumourname[sampleidx], "_subclones.txt", 
  #                                                     sep = ""), rho_psi_file = paste(tumourname[sampleidx], 
  #                                                                                     "_rho_and_psi.txt", sep = ""), gamma_param = platform_gamma)
  # }
  # parallel::stopCluster(clp)
  # if (nsamples > 1) {
  #   print("Assessing mirrored subclonal allelic imbalance (MSAI)")
  #   call_multisample_MSAI(rdsprefix = multisamplehaplotypeprefix, 
  #                         subclonesfiles = paste0(tumourname, "_subclones.txt"), 
  #                         chrom_names = chrom_names, tumournames = tumourname, 
  #                         plotting = T)
  # }
}

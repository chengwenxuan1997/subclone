bbat <- function(analysis = "paired", tumourname, normalname, tumour_data_file, 
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
    if (!skip_phasing[sampleidx]) {
      if (!is.na(externalhaplotypefile) && file.exists(externalhaplotypefile)) {
        externalhaplotypeprefix <- paste0(normalname, 
                                          "_external_haplotypes_chr")
        if (any(!file.exists(paste0(externalhaplotypeprefix, 
                                    1:length(chrom_names), ".vcf")))) {
          print(paste0("Splitting external phasing data from ", 
                       externalhaplotypefile))
          split_input_haplotypes(chrom_names = chrom_names, 
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
          run_haplotyping(chrom = chrom, tumourname = tumourname[sampleidx], 
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
      combine.baf.files(inputfile.prefix = paste(tumourname[sampleidx], 
                                                 "_chr", sep = ""), inputfile.postfix = "_heterozygousMutBAFs_haplotyped.txt", 
                        outputfile = paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", 
                                           sep = ""), chr_names = chrom_names)
    }
    segment.baf.phased(samplename = tumourname[sampleidx], 
                       inputfile = paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", 
                                         sep = ""), outputfile = paste(tumourname[sampleidx], 
                                                                       ".BAFsegmented.txt", sep = ""), prior_breakpoints_file = prior_breakpoints_file, 
                       gamma = segmentation_gamma, phasegamma = phasing_gamma, 
                       kmin = segmentation_kmin, phasekmin = phasing_kmin, 
                       calc_seg_baf_option = calc_seg_baf_option)
    if (nsamples > 1 | write_battenberg_phasing) {
      write_battenberg_phasing(tumourname = tumourname[sampleidx], 
                               SNPfiles = paste0(tumourname[sampleidx], "_alleleFrequencies_chr", 
                                                 chrom_names, ".txt"), imputedHaplotypeFiles = paste0(tumourname[sampleidx], 
                                                                                                      "_impute_output_chr", chrom_names, "_allHaplotypeInfo.txt"), 
                               bafsegmented_file = paste0(tumourname[sampleidx], 
                                                          ".BAFsegmented.txt"), outprefix = paste0(tumourname[sampleidx], 
                                                                                                   "_Battenberg_phased_chr"), chrom_names = chrom_names, 
                               include_homozygous = F)
    }
  }
}

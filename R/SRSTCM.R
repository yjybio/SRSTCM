#' SRSTCM main_func.
#'
#' @param data raw data matrix; if the data is scRNA-seq, genes in rows; cell names in columns; if the data is bulk, genes in rows; sample names in columns.
#' @param cancer short for cancer: "LUNG";"LIHC";"BRCA";"OVCA";"CRC";SKCM".
#' @param type if the data are scRNA-seq, put "sc"; if the data are bulk, put "bulk".
#' @param survival survival data; if the survival data is not empty, survival data is a data frame; sample names, OS, and OS.time in columns.
#' @param multi.sample if scRNA-seq is a single sample, put "FALSE"; if scRNA-seq is multi-sample integration, put "TRUE".
#' @param minGene the minimum number of genes in a cell.
#' @param maxGene the maximum number of genes in a cell.
#' @param maxUMI the maximum UMI count in a cell.
#' @param pctMT mitochondrial gene ratio.
#' @param id.type gene id type: Symbol or Ensemble.
#' @param cell.line if the data are from pure cell line,put "yes"; if cell line data are a mixture of tumor and normal cells, still put "no".
#' @param ngene.chr minimal number of genes per chromosome for cell filtering.
#' @param win.size minimal window sizes for segmentation.
#' @param KS.cut segmentation parameters, input 0 to 1; larger looser criteria.
#' @param distance distance methods include euclidean, and correlation converted distance include pearson and spearman.
#' @param n.cores number of cores for parallel computing.
#'
#' @return if data is scRNA-seq: 1) scRNA-seq quality test results; 2) scRNA-seq cluster map; 3) aneuploid/diploid prediction results; 4)senescent tumor cells prediction results; 5)annotated cluster diagram of senescent tumor cells; 6)differential genes in senescent and non-senescent tumor cell clusters; 7)expression of the top 10 differentially expressed genes in senescent and non-senescent tumor cell clusters; if data is bulk: 1)the senescence score results; 2)survival outcomes in the high and low senescence score groups.
#' @export
#'
#' @examples
#' test.SRSTCM <- SRSTCM(data = data, cancer = "LUNG", type = "sc", survival = survival, multi.sample = "FALSE", minGene = 500, maxGene = 5000, maxUMI = 40000, pctMT = 20, id.type = "S", cell.line = "no", ngene.chr = 5,  win.size = 25, KS.cut = 0.15, distance = "euclidean", n.cores = 1)
SRSTCM <- function (data, cancer, type, survival, multi.sample = "FALSE",
                    minGene = 500, maxGene = 5000, maxUMI = 40000,
                    pctMT = 20, id.type = "S", cell.line = "no", ngene.chr = 5,
                    win.size = 25, KS.cut = 0.15, distance = "euclidean", n.cores = 1)
{
  CRC.sen.gene <- c("ANG", "ANGPTL4", "AREG", "AdataL", "BMP6", "C3", "CCL20", "CCL26", "CD55", "CD9", "CSF1", "CSF2", "CTSB", "CdataCL1", "CdataCL16", "CdataCL2", "CdataCL3", "CdataCL8", "DKK1", "EDN1", "EGFR", "EREG", "ETS2", "FAS", "FGF1", "FGF7", "GDF15", "GEM", "GMFG", "HMGB1", "ICAM1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL15", "IL18", "IL1B", "IL6ST", "INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "LCP1", "MIF", "MMP1", "MMP14", "NAP1L4", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGES", "RPS6KA5", "SCAMP4", "SEMA3F", "SERPINE1", "SERPINE2", "SPP1", "SPdata", "TIMP2", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VGF", "WNT16")
  CRC.sen.average <- c(1.525989,0.1722437,3.803822,4.359717,0.4387834,0.7239708,2.75195,1.2432,1.462746,8.47853,0.3163165,0.4596274,7.202994,2.181346,2.475049,0.5578473,2.370216,3.995862,3.505363,0.786515,0.4442448,1.986943,5.117946,0.8846947,0.1656492,0.01382584,7.433819,0.3222421,3.931369,8.69461,2.052468,2.846775,4.944452,2.49919,7.261223,0.07010185,6.790036,0.1066997,1.30141,0.08425097,1.862054,0.06194288,2.433603,1.223618,0.7712733,2.921095,3.961382,5.615334,6.965025,9.196636,0.0980821,3.288823,6.184831,0.1279435,1.257418,3.205868,1.643292,3.623119,3.660525,8.240554,0.8995773,2.280486,5.269145,5.271834,4.550359,5.071742,0.2233027,1.273417,6.372992,0.2479643,1.050233,5.418115,4.801694,4.130594,4.497571,1.895459,0.01087502)
  CRC.no.sen.average <- c(1.762392,0.05609277,1.655296,3.832673,0.2042809,0.06887844,2.550213,0.1927539,0.3610346,8.517132,0.1285929,0.05015286,6.437942,1.995801,1.276756,0.4733429,2.141877,3.210712,2.105524,1.89567,0.0335546,0.5985007,5.467061,0.505045,0.1136677,0.01847663,6.786217,0.02735736,3.769336,9.462114,1.156614,1.196588,5.402516,1.057567,7.410258,0.1173942,6.115179,0.09703373,1.109096,0.02790151,1.555718,0.01947779,2.176475,0.8567837,0.1280684,2.024906,2.681027,5.799405,6.901019,9.756015,0.0322134,2.300517,6.443458,0.1359118,0.6191497,3.515789,0.8038596,2.240058,2.320443,8.555899,0.5007245,2.076116,5.102469,5.035261,2.856537,3.791328,0.1801424,0.5049134,5.829992,0.1471127,1.236701,5.248839,4.825674,4.298734,4.395311,1.285141,0.02260939)
  LUNG.sen.gene <- c("CCL5","FGF1","ANGPT1","RPS6KA5","FGF7","CD55","IGFBP5","NAP1L4","AREG","PECAM1","PAPPA","EDN1","EGFR","MMP2","SPP1","SERPINE1","CXCL2","EREG","MMP14","HMGB1","CTSB","CD9","IGFBP7","ETS2","JUN","IGFBP4","PLAT","TUBGCP2","ICAM1","AXL","IGFBP2","CXCL8","TIMP2","IQGAP2","TNFRSF1B","CXCL12","IGFBP6","MMP9","CCL4","GMFG","FGF2","CXCL1","GEM","MMP1","CXCL10","MMP12","DKK1","FAS","IL6ST","TNFRSF11B","ITGA2","IL1B","PIGF","ANG","CSF2RB","IL6","BMP2","IGFBP1","CCL20","PLAU","VGF","WNT2","MMP10","MMP3","ITPKA","IL15","BMP6","TNFRSF10C","NRG1","EGF","IL18","PTGER2","IL7","SEMA3F","CST4","KITLG","CSF1","TNF","PTGES","IL10","TNFRSF1A","IL2","CXCL3","IL1A","ESM1","LCP1","PGF","VEGFC","HGF","IGFBP3","INHA","CSF2","VEGFA","PLAUR","PTBP1","SERPINB4","SERPINE2","SCAMP4","CCL8","CCL2","C3","MIF","ANGPTL4","SPX","WNT16","GDF15","CXCL16","CCL26")
  LUNG.sen.average <- c(4.516653,7.11344,7.028329,6.440879,4.990518,10.17978,11.39231,9.885077,11.50175,6.0829,12.3509,10.12917,9.21291,11.45692,6.356315,12.44274,12.56048,13.22759,12.21906,13.06438,13.48016,12.21544,13.1311,10.63429,11.84318,9.818679,13.1952,8.508703,11.29588,10.33655,3.869751,15.54269,12.71481,6.544546,9.854555,5.470181,10.00026,8.106235,5.574034,7.032608,12.58717,15.22726,11.47288,15.79085,11.31467,11.11047,9.771289,10.67309,11.54399,11.41186,13.87191,15.69128,10.78496,6.565153,2.922089,15.15221,13.18889,6.422209,12.91271,15.02163,6.694204,3.534783,10.25335,15.65018,7.234016,7.67153,4.025343,5.504512,10.83957,5.348262,3.090463,7.28894,4.594162,9.250072,10.13507,10.79707,6.281226,3.942871,6.95573,4.465434,10.76052,2.58193,14.45098,12.35661,9.901724,5.267071,8.409468,12.80502,10.54311,10.55381,6.124658,12.31862,11.42296,12.86923,10.75909,10.99261,12.5205,7.319861,4.482969,10.7125,9.683254,14.66438,9.486766,5.996619,8.344061,12.86749,3.008246,5.124981)
  LUNG.no.sen.average <- c(3.499048,8.224109,9.342742,5.690249,7.162573,10.28033,12.62224,9.679599,7.278845,5.246589,9.881582,8.70109,9.247568,11.09202,4.714892,12.18619,4.469914,8.457313,10.85492,13.67238,13.15685,12.0203,14.23808,10.48947,11.57702,10.43796,12.71669,8.971236,8.446585,12.7647,4.745125,9.144794,13.50005,4.894808,8.10599,8.407087,12.12614,7.179001,2.349387,6.006501,10.43933,6.589595,10.39668,14.20804,7.829934,5.153116,12.25673,10.7071,11.63145,12.16819,12.95435,7.317924,11.25958,7.512106,3.11633,10.73093,10.29341,4.356221,2.207406,14.04162,5.794648,7.727369,6.665612,10.5796,4.954662,9.211342,4.828784,5.01651,7.019959,5.462573,4.926468,10.86016,3.015609,9.19806,6.684621,10.33513,7.722617,5.423998,6.509167,2.11633,11.31166,1.938872,5.45999,5.183628,6.353673,4.607062,8.568338,12.1474,11.17033,11.8558,5.113868,5.965854,9.220687,11.31227,10.57356,1.572023,12.49333,7.357542,3.669925,10.72751,4.686955,14.23449,5.327622,3.720953,9.196812,9.705761,3.648687,5.91517)
  SKCM.sen.gene <- c("AREG","ITGA2","PLAUR","MMP12","IGFBP2","MMP1","SEMA3F","PLAU","NAP1L4","HMGB1","PTGES","EDN1","IL1A","IQGAP2","ANGPT1","PTBP1","BMP6","IL6","MMP9","GMFG","GEM","SERPINE1","MMP10","HGF","CSF2","TNFRSF10C","IL1B","TNFRSF1A","GDF15","VEGFA","EGF","SERPINE2","IGFBP4","LCP1","ESM1","AXL","CXCL1","IGFBP3","CXCL2","KITLG","ANG","FGF7","CXCL3","CD9","TIMP2","FAS","DKK1","FGF1","PIGF","SCAMP4","TNFRSF11B","PGF","CD55","IGFBP5","JUN","CSF1","PLAT","IL6ST","CCL2","MMP3","EREG","RPS6KA5","SPP1","ITPKA","PTGER2","ETS2","PECAM1","CCL20","IL15","CXCL8","VEGFC","CTSB","EGFR","NRG1","WNT16","FGF2","IL18","CXCL10","IL7","TUBGCP2","CXCL16","CCL5","TNFRSF1B","MMP14","MMP2","ANGPTL4","VGF","IGFBP6","C3","IL10","CSF2RB","PAPPA","BMP2","INHA","IGFBP7","MIF","ICAM1")
  SKCM.sen.average <- c(4.534673,4.502458,5.453566,0.9127044,5.208838,2.623085,1.314262,0.1447677,5.044333,6.7535,2.453942,0.1133612,4.195601,0.5982176,2.43401,6.65909,2.215571,0.1392059,0.7155627,1.177499,3.722463,1.050253,1.472719,0.8312834,1.070051,1.664389,2.80752,6.671775,3.13186,4.315163,0.7078222,11.11836,5.936078,0.3269915,1.130398,4.841613,6.905172,6.913744,3.284899,0.7505059,1.195812,1.16542,3.793877,7.856049,5.724522,4.142248,5.091707,2.066186,4.590378,4.313709,2.354252,1.153786,7.906183,3.987674,5.239003,2.889257,6.467422,5.05848,6.554183,0.8997989,0.07775439,1.412717,3.108323,0.2801438,0.104704,1.785376,0.2904251,0.8906483,0.6419515,6.958457,2.03181,6.438099,1.141349,0.7936378,0.5161879,2.121318,0.4582779,1.617637,1.317504,4.928349,3.340874,1.437157,0.6806705,6.087502,7.148653,2.321936,6.098573,4.544428,6.851433,1.207128,0.02364627,2.149832,2.785605,1.919894,5.460785,8.596628,6.882476)
  SKCM.no.sen.average <- c(4.853115,4.88515,6.676557,0.1607168,5.690992,3.875924,2.140709,0.2115144,5.181754,7.854619,4.933307,0.140154,6.249707,2.025673,2.344408,7.343832,3.257839,0.5528465,4.40633,0.8628525,4.473123,1.830588,2.902651,0.5588978,5.294609,1.069686,4.379294,6.822238,2.764869,5.000979,0.3278449,11.82118,5.148068,0.2613215,4.611981,5.697242,8.813927,9.314631,5.662279,0.3652238,2.593964,1.026008,6.568657,7.354802,5.503427,3.971878,4.680936,0.5091236,3.815126,4.051171,2.94625,1.110921,7.940614,4.858929,4.082465,2.144437,6.024994,5.183368,3.993266,0.08187245,0.126901,2.185619,2.459415,1.597886,0.1106601,2.717263,0.8306982,2.63477,1.222168,9.117962,3.417052,6.486434,2.494439,0.7506097,0.6296158,2.056509,0.543558,0.9243875,1.659288,5.014145,3.084058,1.186863,1.713507,6.197054,6.440989,3.818006,7.663923,6.207481,6.829956,1.927612,0.01899988,1.620843,3.744709,3.241613,6.090823,9.355998,7.199227)
  LIHC.sen.gene <- c("CCL5","FGF1","ANGPT1","RPS6KA5","FGF7","CD55","IGFBP5","NAP1L4","AREG","PECAM1","PAPPA","EDN1","EGFR","MMP2","SPP1","SERPINE1","CXCL2","EREG","MMP14","HMGB1","CTSB","CD9","IGFBP7","ETS2","JUN","IGFBP4","PLAT","TUBGCP2","ICAM1","AXL","IGFBP2","CXCL8","TIMP2","IQGAP2","TNFRSF1B","CXCL12","IGFBP6","MMP9","CCL4","GMFG","FGF2","CXCL1","GEM","MMP1","CXCL10","MMP12","DKK1","FAS","IL6ST","TNFRSF11B","ITGA2","IL1B","PIGF","ANG","CSF2RB","IL6","BMP2","IGFBP1","CCL20","PLAU","VGF","WNT2","MMP10","MMP3","ITPKA","IL15","BMP6","TNFRSF10C","NRG1","EGF","IL18","PTGER2","IL7","SEMA3F","CST4","KITLG","CSF1","TNF","PTGES","IL10","TNFRSF1A","IL2","CXCL3","IL1A","ESM1","LCP1","PGF","VEGFC","HGF","IGFBP3","INHA","CSF2","VEGFA","PLAUR","PTBP1","SERPINB4","SERPINE2","SCAMP4","CCL8","CCL2","C3","MIF","ANGPTL4","SPX","WNT16","GDF15","CXCL16","CCL26")
  LIHC.sen.average <- c(3.015722,2.996278,4.192444,3.138167,2.353778,7.925167,3.403278,5.136767,5.012,4.6911,3.213591,7.075167,7.164,3.041611,3.581167,4.865722,4.911167,3.078833,4.923083,9.5275,7.783433,5.148611,3.415333,5.808444,6.560417,4.125667,2.837667,6.118944,7.463444,3.399583,4.143667,9.445,6.320333,5.969833,3.19,3.122833,5.862833,3.902667,3.4005,3.373833,2.949083,9.196167,5.716833,2.7625,5.763167,2.6435,7.904833,2.779542,5.757262,6.63075,6.61125,3.7965,7.33475,8.134667,2.549167,3.261,10.33833,3.76525,12.35367,4.898417,4.716167,3.0565,2.842667,2.580833,6.526833,2.947083,3.442167,3.627292,3.511611,5.901833,6.944333,3.032667,3.820167,5.132778,3.709333,6.661056,5.463125,3.475,4.503,2.515167,7.948167,2.324667,7.075333,3.02375,2.473167,6.343,7.1445,2.694333,3.007767,8.917167,3.196333,3.97325,7.375667,4.343889,5.782292,2.412333,6.058083,5.098333,2.5575,2.707667,12.05383,11.06767,4.382083,6.840417,3.0795,5.041278,7.094333,2.9715)
  LIHC.no.sen.average <- c(3.049556,3.043944,4.585722,2.909944,2.411778,6.837889,3.447194,5.3787,4.037722,4.8054,3.247652,5.666333,7.007056,3.051444,3.040167,5.427722,4.357722,3.076,4.895792,10.34823,7.4369,4.503889,3.454611,5.761167,6.194208,3.956833,2.729,6.406222,7.012,3.3895,3.6775,7.876583,4.384833,6.0465,3.092167,3.301917,5.143333,3.869167,3.218167,3.421,3.382167,8.581,5.453,2.5355,4.535667,2.611833,9.956833,2.829917,5.876238,7.354667,5.801833,3.675833,7.076583,7.317333,2.795167,3.065833,10.11883,3.249083,11.61117,5.001917,4.834167,3.101667,2.666,2.743833,7.2255,2.89825,3.450889,3.66725,3.583222,3.715167,6.464167,2.571333,4.253167,5.121,3.791167,6.930944,4.760375,3.495167,4.39625,2.426833,7.3685,2.522333,5.326833,2.945917,2.505333,5.714167,6.7295,2.755667,3.057033,8.117167,3.343667,4.070583,7.16825,4.220056,5.990667,2.4525,5.927917,5.812167,2.576667,2.733833,11.68117,11.30117,4.515667,5.1965,3.130417,5.044333,6.260667,2.915667)
  BRCA.sen.gene <- c("ANG","ANGPTL4","AREG","AXL","BMP6","C3","CCL26","CCL5","CD55","CD9","CSF1","CTSB","CXCL12","CXCL16","DKK1","EDN1","EGF","EGFR","EREG","ETS2","FAS","GDF15","GEM","GMFG","HMGB1","ICAM1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IL15","IL18","IL1A","IL1B","IL6ST","IL7","INHA","ITGA2","ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP14","MMP9","NAP1L4","PAPPA","PGF","PIGF","PLAT","PLAU","PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SEMA3F","SERPINE1","SERPINE2","TIMP2","TNFRSF10C","TNFRSF11B","TNFRSF1A","TUBGCP2","VEGFA","VEGFC","VGF","WNT16")
  BRCA.sen.average <- c(1.633935,3.287,6.196108,0.02446049,0.3519869,4.322028,1.288669,0.06309963,4.509007,5.370831,1.260364,4.868937,1.462404,3.486612,6.050021,2.158759,0.5092545,1.583384,0.3695031,4.207032,1.536849,6.507353,0.4739225,0.1090695,3.821826,3.607362,4.496967,3.657419,5.949303,2.449754,0.8294141,0.04602155,1.411269,0.02621634,0.02674078,2.948694,0.07225448,2.164955,3.823793,1.623743,4.061699,4.144362,2.338374,0.7440392,0.8082699,0.3017092,0.728538,4.613928,0.03118674,0.1264931,1.852731,0.1120781,0.6917155,1.209772,6.411088,1.812422,5.616917,0.453628,3.128922,3.991994,0.04694525,1.32766,3.626739,2.826806,2.665809,6.257189,3.619016,2.466752,1.223284,3.788017,0.04263861)
  BRCA.no.sen.average <- c(0.7291369,1.363872,5.552337,0.03305695,0.2889156,0.1363496,0.2429305,0.04607023,3.921063,6.744796,0.1960821,4.558778,4.600158,2.88801,6.466418,1.551508,0.2190206,1.235165,0.2299785,2.36754,1.023388,3.862153,0.4560649,0.07461124,4.223604,0.9037696,4.132839,2.0114,7.790752,5.367937,1.514998,0.02395574,0.4442801,0.1237429,0.02714217,2.437825,0.06614268,1.200897,2.808807,1.606503,3.276966,3.750201,1.782997,0.6667997,0.06471584,0.361439,0.8862556,4.790388,0.02453238,0.2592028,1.829622,0.004207816,0.2188561,0.8193248,6.460266,1.138417,1.370382,0.623867,3.408061,4.627124,0.2632277,0.3873274,3.423886,1.674419,2.030284,5.206165,3.759438,1.855281,1.036696,3.103537,0.05101531)
  OVCA.sen.gene <- c("CXCL8","CXCL1","PLAU","CXCL2","CD55","PLAT","IGFBP2","C3","JUN","PLAUR","SPP1","PTBP1","SERPINE2","CD9","ITGA2","IGFBP7","AREG","SEMA3F","PTGES","CTSB","HMGB1","CXCL3","NRG1","ETS2","EDN1","IL6","IGFBP3","MMP10","VEGFA","IGFBP6","IL18","IGFBP5","ICAM1","ANGPT1","CCL20","FGF2","CXCL10","IGFBP4","DKK1","EGF","IL15","CSF1","IL6ST","PTGER2","MMP9","TIMP2","LCP1","SPX","RPS6KA5","MMP2","IQGAP2","GMFG","SERPINE1","EREG","TNFRSF1A","KITLG","CXCL16","EGFR","ITPKA","PGF","SCAMP4","PIGF","TUBGCP2","FAS","GDF15")
  OVCA.sen.average <- c(2.438247,7.535038,4.541467,1.14888,5.035944,1.753331,2.721673,5.896031,4.655341,2.267511,1.877889,4.701353,0.730662,7.199818,0.4009275,2.84992,4.493424,2.184036,2.387924,7.876397,3.981276,2.034359,1.331385,5.011053,1.433098,0.4091184,4.279176,0.8839654,1.856869,1.095935,2.969512,1.808716,1.321868,0.4157197,2.980465,0.2443182,1.008442,5.158154,0.2270614,0.09455035,0.4837918,0.6258139,3.113056,0.465292,1.965323,5.895799,1.435189,2.135363,0.335606,0.8764409,0.1255488,1.105786,0.4392672,0.298693,4.176954,0.3236572,4.110859,2.656538,0.5300761,0.3197893,2.046181,1.223216,3.63654,2.821648,2.640621)
  OVCA.no.sen.average <- c(5.530194,10.03251,6.684943,4.678461,3.569932,3.042536,4.709664,4.587779,5.765207,3.578608,0,5.51334,1.666505,8.218488,1.178418,4.231154,3.396588,3.272751,3.667746,8.710727,4.724384,2.975675,1.9341,4.456846,2.273826,1.440952,3.871161,0.3258267,1.464558,1.862385,3.427362,1.440536,1.712209,0.1510116,3.542978,0.1078651,0.2815664,5.505053,0.5445676,0.2461861,0.6756182,1.013597,2.862678,0.8625667,2.339958,6.052229,1.655585,1.901131,0.2106882,0.9312557,0.09386151,1.755048,0.4278183,0.2372753,4.254538,0.3370817,3.920511,2.706014,0.4933547,0.4629519,2.092599,1.109999,3.715668,2.72128,2.773165)
  q <- function (x) {
    (x + 1) / ( sum(x) + length(x) )
  }
  argmax <- function (p, q) {
    result <- 0
    for (i in 1 : length(p)) {
      result <- result + ( q[i] * log(p[i]) - q[i] * log(q[i]) )
    }
    return(result)
  }
  if (type == "sc") {
    print("Step1: Single cell RNA-seq data reading and quality control ...")
    if (multi.sample == "TRUE") {
      data <- subset(data, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene &
                       nFeature_RNA < maxGene & percent.mt < pctMT)
      data <- NormalizeData(data) %>% FindVariableFeatures(nfeatures =
                                                             ceiling(nrow(data@assays$RNA) * 0.3)) %>% ScaleData()
      data <- RunPCA(data, features = VariableFeatures(object = data), reduction.name = "pca")
      data <- RunUMAP(data, reduction = "pca", dims = 1:10, reduction.name = "umap_naive")
      data <- RunHarmony(data, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
      data <- RunUMAP(data, reduction = "harmony", dims = 1:30,reduction.name = "umap")
      data <- FindNeighbors(data, reduction = "harmony", dims = 1:30)
      data <- FindClusters(data, res = 0.5)
    }
    else {
      data <- subset(data, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene &
                       nFeature_RNA < maxGene & percent.mt < pctMT)
      data <- NormalizeData(data) %>% FindVariableFeatures(nfeatures =
                                                             ceiling(nrow(data@assays$RNA) * 0.3)) %>% ScaleData()
      data <- RunPCA(data, features = VariableFeatures(object = data), reduction.name = "pca")
      data <- RunUMAP(data, reduction = "pca", dims = 1:10, reduction.name = "umap")
      data <- FindNeighbors(data, dims = 1:30)
      data <- FindClusters(data, res = 0.5)
    }
    saveRDS(data, file = "quality_testing_result.rds")
    pdf(paste(cancer, "_seurat_clusters.pdf", sep = ""), height = 6, width = 8)
    p1 <- DimPlot(data, reduction = "umap", label = TRUE, group.by = "seurat_clusters", label.size = 3)
    print(p1)
    dev.off()
    print("Step2: Predictive tumor cells ...")
    tumor.cell <- as.matrix(data@assays$RNA@counts)
    copykat.test <- copykat(rawmat = tumor.cell, sam.name = cancer)
    sample.name <- paste(cancer, "_copykat_", sep = "")
    res <- paste(sample.name, "prediction.txt", sep = "")
    copykat.prediction <- read.table(res)
    copykat.prediction <- copykat.prediction[which(copykat.prediction[, 2] == "aneuploid"), ]
    copykat.prediction <- data@assays$RNA@counts[ ,which(colnames(data@assays$RNA@counts) %in% copykat.prediction[, 1])]
    print("Step3: Final prediction ...")
    if (cancer == "CRC") {
      copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% CRC.sen.gene), ]
      dir <- which(CRC.sen.gene %in% rownames(copykat.prediction))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        copykat.prediction <- copykat.prediction[match(CRC.sen.gene[dir], rownames(copykat.prediction)), ]
        CRC.sen.average <- CRC.sen.average[dir]
        CRC.no.sen.average <- CRC.no.sen.average[dir]
        p.no.sen <- ( CRC.no.sen.average + 1 ) / ( sum(CRC.no.sen.average) + length(CRC.no.sen.average) )
        p.sen <- ( CRC.sen.average + 1 ) / ( sum(CRC.sen.average) + length(CRC.sen.average) )
        q.CRC <- apply(copykat.prediction, 2, q)
        senescence <- c()
        for (i in 1 : ncol(q.CRC)) {
          if (argmax(p.sen, q.CRC[, i]) > argmax(p.no.sen, q.CRC[, i])) {
            senescence <- c(senescence, colnames(q.CRC)[i])
          }
        }
        write.table(senescence, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "LUNG") {
      copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% LUNG.sen.gene), ]
      dir <- which(LUNG.sen.gene %in% rownames(copykat.prediction))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        copykat.prediction <- copykat.prediction[match(LUNG.sen.gene[dir], rownames(copykat.prediction)), ]
        LUNG.sen.average <- LUNG.sen.average[dir]
        LUNG.no.sen.average <- LUNG.no.sen.average[dir]
        p.no.sen <- ( LUNG.no.sen.average + 1 ) / ( sum(LUNG.no.sen.average) + length(LUNG.no.sen.average) )
        p.sen <- ( LUNG.sen.average + 1 ) / ( sum(LUNG.sen.average) + length(LUNG.sen.average) )
        q.LUNG <- apply(copykat.prediction, 2, q)
        senescence <- c()
        for (i in 1 : ncol(q.LUNG)) {
          if (argmax(p.sen, q.LUNG[, i]) > argmax(p.no.sen, q.LUNG[, i])) {
            senescence <- c(LUNG.sen.average, colnames(q.LUNG)[i])
          }
        }
        write.table(senescence, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "SKCM") {
      copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% SKCM.sen.gene), ]
      dir <- which(SKCM.sen.gene %in% rownames(copykat.prediction))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        copykat.prediction <- copykat.prediction[match(SKCM.sen.gene[dir], rownames(copykat.prediction)), ]
        SKCM.sen.average <- SKCM.sen.average[dir]
        SKCM.no.sen.average <- SKCM.no.sen.average[dir]
        p.no.sen <- ( SKCM.no.sen.average + 1 ) / ( sum(SKCM.no.sen.average) + length(SKCM.no.sen.average) )
        p.sen <- ( SKCM.sen.average + 1 ) / ( sum(SKCM.sen.average) + length(SKCM.sen.average) )
        q.SKCM <- apply(copykat.prediction, 2, q)
        senescence <- c()
        for (i in 1 : ncol(q.SKCM)) {
          if (argmax(p.sen, q.SKCM[, i]) > argmax(p.no.sen, q.SKCM[, i])) {
            senescence <- c(senescence, colnames(q.SKCM)[i])
          }
        }
        write.table(senescence, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "LIHC") {
      copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% LIHC.sen.gene), ]
      dir <- which(LIHC.sen.gene %in% rownames(copykat.prediction))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) >1) {
        copykat.prediction <- copykat.prediction[match(LIHC.sen.gene[dir], rownames(copykat.prediction)), ]
        LIHC.sen.average <- LIHC.sen.average[dir]
        LUNG.no.sen.average <- LUNG.no.sen.average[dir]
        p.no.sen <- ( LUNG.no.sen.average + 1 ) / ( sum(LUNG.no.sen.average) + length(LUNG.no.sen.average) )
        p.sen <- ( LIHC.sen.average + 1 ) / ( sum(LIHC.sen.average) + length(LIHC.sen.average) )
        q.LIHC <- apply(copykat.prediction, 2, q)
        senescence <- c()
        for (i in 1 : ncol(q.LIHC)) {
          if (argmax(p.sen, q.LIHC[, i]) > argmax(p.no.sen, q.LIHC[, i])){
            senescence <- c(senescence, colnames(q.LIHC)[i])
          }
        }
        write.table(senescence, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "BRCA") {
      copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% BRCA.sen.gene), ]
      dir <- which(BRCA.sen.gene %in% rownames(copykat.prediction))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        copykat.prediction <- copykat.prediction[match(BRCA.sen.gene[dir], rownames(copykat.prediction)), ]
        BRCA.sen.average <- BRCA.sen.average[dir]
        BRCA.no.sen.average <- BRCA.no.sen.average[dir]
        p.no.sen <- ( BRCA.no.sen.average + 1 ) / ( sum(BRCA.no.sen.average) + length(BRCA.no.sen.average) )
        p.sen <- ( BRCA.sen.average + 1 ) / ( sum(BRCA.sen.average) + length(BRCA.sen.average) )
        q.BRCA <- apply(copykat.prediction, 2, q)
        senescence <- c()
        for (i in 1 : ncol(q.BRCA)) {
          if (argmax(p.sen, q.BRCA[, i]) > argmax(p.no.sen, q.BRCA[, i])) {
            senescence <- c(senescence, colnames(q.BRCA)[i])
          }
        }
        write.table(senescence, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "OVCA") {
      copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% OVCA.sen.gene), ]
      dir <- which(OVCA.sen.gene %in% rownames(copykat.prediction))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        copykat.prediction <- copykat.prediction[match(OVCA.sen.gene[dir], rownames(copykat.prediction)), ]
        OVCA.sen.average <- OVCA.sen.average[dir]
        OVCA.no.sen.average <- OVCA.no.sen.average[dir]
        p.no.sen <- ( OVCA.no.sen.average + 1 ) / ( sum(OVCA.no.sen.average) + length(OVCA.no.sen.average) )
        p.sen <- ( OVCA.sen.average + 1 ) / ( sum(OVCA.sen.average) + length(OVCA.sen.average) )
        q.OVCA <- apply(copykat.prediction, 2, q)
        senescence <- c()
        for (i in 1 : ncol(q.OVCA)) {
          if (argmax(p.sen, q.OVCA[, i]) > argmax(p.no.sen, q.OVCA[, i])) {
            senescence <- c(senescence, colnames(q.OVCA)[i])
          }
        }
        write.table(senescence, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
      }
    }
    data@meta.data$cell.type <- data@meta.data$seurat_clusters
    data@meta.data$cell.type <- as.character(data@meta.data$cell.type)
    data@meta.data$cell.type[which(rownames(data@meta.data) %in% senescence)] <- "STCs"
    data@meta.data$cell.type <- as.factor(data@meta.data$cell.type)
    pdf(paste(cancer, "_note_STCs.pdf", sep = ""), height = 6, width = 8)
    p2 <- DimPlot(data, reduction = "umap", label = TRUE, group.by = "cell.type", label.size = 3)
    print(p2)
    dev.off()
    data@meta.data$cell.type <- as.character(data@meta.data$cell.type)
    no.senescence <- colnames(copykat.prediction)[-which(colnames(copykat.prediction) %in% senescence)]
    data@meta.data$cell.type[which(rownames(data@meta.data) %in% no.sen)] <- "NSTCs"
    data@meta.data$cell.type <- as.factor(data@meta.data$cell.type)
    Idents(data) <- data@meta.data$cell.type
    mydeg <- FindMarkers(data, ident.1 = 'STCs', ident.2 = 'NSTCs', verbose = FALSE, test.use = 'wilcox', min.pct = 0.1)
    write.csv(mydeg, paste(cancer, "_STCs_NSTCs_deg.csv", sep = ""), row.names = TRUE)
    top10 <- mydeg %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
    pdf(paste(cancer, "_top10__STCs_dge.pdf", sep = ""), height = 6, width = 8)
    p3 <- VlnPlot(data, features = top10, split.by = 'cell.type', idents = 'STCs')
    print(p3)
    dev.off()
    pdf(paste(cancer, "_top10__NSTCs_dge.pdf", sep = ""), height = 6, width = 8)
    p4 <- VlnPlot(data, features = top10, split.by = 'cell.type', idents = 'NSTCs')
    print(p4)
    dev.off()
  }
  if (type == "bulk") {
    print("Read and analyze bulk RNA-seq Data ...")
    name <- rownames(data)
    data <- apply(data, 2, as.numeric)
    rownames(data) <- name
    if (cancer == "CRC") {
      data <- data[which(rownames(data) %in% CRC.sen.gene), ]
      dir <- which(CRC.sen.gene %in% rownames(data))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        data <- data[match(CRC.sen.gene[dir], rownames(data)), ]
        CRC.sen.average <- CRC.sen.average[dir]
        p.sen <- ( CRC.sen.average + 1 ) / ( sum(CRC.sen.average) + length(CRC.sen.average) )
        q.CRC <- apply(data, 2, q)
        Senescence.score <- matrix(nrow = ncol(q.CRC), ncol = 2)
        for (i in 1 : ncol(q.CRC)) {
          Senescence.score[i, 1] <- colnames(q.CRC)[i]
          Senescence.score[i, 2] <- argmax(p.sen, q.CRC[, i])
        }
        colnames(Senescence.score) <- c("sample", "senescence_score")
        Senescence.score <- as.data.frame(Senescence.score)
        write.table(Senescence.score, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "LUNG") {
      data <- data[which(rownames(data) %in% LUNG.sen.gene), ]
      dir <- which(LUNG.sen.gene %in% rownames(data))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        data <- data[match(LUNG.sen.gene[dir], rownames(data)), ]
        LUNG.sen.average <- LUNG.sen.average[dir]
        p.sen <- ( LUNG.sen.average + 1 ) / ( sum(LUNG.sen.average) + length(LUNG.sen.average) )
        q.LUNG <- apply(data, 2, q)
        Senescence.score <- matrix(nrow = ncol(q.LUNG), ncol = 2)
        for (i in 1 : ncol(q.LUNG)) {
          Senescence.score[i, 1] <- colnames(q.LUNG)[i]
          Senescence.score[i, 2] <- argmax(p.sen, q.LUNG[, i])
        }
        colnames(Senescence.score) <- c("sample", "senescence_score")
        Senescence.score <- as.data.frame(Senescence.score)
        write.table(Senescence.score, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "SKCM") {
      data <- data[which(rownames(data) %in% SKCM.sen.gene), ]
      dir <- which(SKCM.sen.gene %in% rownames(data))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        data <- data[match(SKCM.sen.gene[dir], rownames(data)), ]
        SKCM.sen.average <- SKCM.sen.average[dir]
        p.sen <- ( SKCM.sen.average + 1 ) / ( sum(SKCM.sen.average) + length(SKCM.sen.average) )
        q.SKCM <- apply(data, 2, q)
        Senescence.score <- matrix(nrow = ncol(q.SKCM), ncol = 2)
        for (i in 1 : ncol(q.SKCM)) {
          Senescence.score[i, 1] <- colnames(q.SKCM)[i]
          Senescence.score[i, 2] <- argmax(p.sen, q.SKCM[, i])
        }
        colnames(Senescence.score) <- c("sample", "senescence_score")
        Senescence.score <- as.data.frame(Senescence.score)
        write.table(Senescence.score, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "LIHC") {
      data <- data[which(rownames(data) %in% LIHC.sen.gene), ]
      dir <- which(LIHC.sen.gene %in% rownames(data))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) >1) {
        data <- data[match(LIHC.sen.gene[dir], rownames(data)), ]
        LIHC.sen.average <- LIHC.sen.average[dir]
        p.sen <- ( LIHC.sen.average + 1 ) / ( sum(LIHC.sen.average) + length(LIHC.sen.average) )
        q.LIHC <- apply(data, 2, q)
        Senescence.score <- matrix(nrow = ncol(q.LIHC), ncol = 2)
        for (i in 1 : ncol(q.LIHC)) {
          Senescence.score[i, 1] <- colnames(q.LIHC)[i]
          Senescence.score[i, 2] <- argmax(p.sen, q.LIHC[, i])
        }
        colnames(Senescence.score) <- c("sample", "senescence_score")
        Senescence.score <- as.data.frame(Senescence.score)
        write.table(Senescence.score, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "BRCA") {
      data <- data[which(rownames(data) %in% BRCA.sen.gene), ]
      dir <- which(BRCA.sen.gene %in% rownames(data))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        data <- data[match(BRCA.sen.gene[dir], rownames(data)), ]
        BRCA.sen.average <- BRCA.sen.average[dir]
        p.sen <- ( BRCA.sen.average + 1 ) / ( sum(BRCA.sen.average) + length(BRCA.sen.average) )
        q.BRCA <- apply(data, 2, q)
        Senescence.score <- matrix(nrow = ncol(q.BRCA), ncol = 2)
        for (i in 1 : ncol(q.BRCA)) {
          Senescence.score[i, 1] <- colnames(q.BRCA)[i]
          Senescence.score[i, 2] <- argmax(p.sen, q.BRCA[, i])
        }
        colnames(Senescence.score) <- c("sample", "senescence_score")
        Senescence.score <- as.data.frame(Senescence.score)
        write.table(Senescence.score, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
    if (cancer == "OVCA") {
      data <- data[which(rownames(data) %in% OVCA.sen.gene), ]
      dir <- which(OVCA.sen.gene %in% rownames(data))
      if (length(dir) == 0 || length(dir) == 1)
        stop("too few senescence-related genes;cannot be calculated")
      if (length(dir) > 1) {
        data <- data[match(OVCA.sen.gene[dir], rownames(data)), ]
        OVCA.sen.average <- OVCA.sen.average[dir]
        p.sen <- ( OVCA.sen.average + 1 ) / ( sum(OVCA.sen.average) + length(OVCA.sen.average) )
        q.OVCA <- apply(data, 2, q)
        Senescence.score <- matrix(nrow = ncol(q.OVCA), ncol = 2)
        for (i in 1 : ncol(q.OVCA)) {
          Senescence.score[i, 1] <- colnames(q.OVCA)[i]
          Senescence.score[i, 2] <- argmax(p.sen, q.OVCA[, i])
        }
        colnames(Senescence.score) <- c("sample", "senescence_score")
        Senescence.score <- as.data.frame(Senescence.score)
        write.table(Senescence.score, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
    if (!is.null(survival)) {
      Senescence.score$senescence_score <- as.numeric(Senescence.score$senescence_score)
      Senescence.score <- Senescence.score[match(sort(Senescence.score$senescence_score), Senescence.score$senescence_score), ]
      survival <- survival[survival$sample %in% Senescence.score$sample, ]
      survival <- survival[match(Senescence.score$sample, survival$sample), ]
      survival$senescence_score <- Senescence.score$senescence_score
      count <- ceiling(nrow(survival) / 2)
      class <- c(rep("low_score", count), rep("high_score", (nrow(survival) - count)))
      survival$class <- class
      fit <- survfit(Surv(OS.time, OS) ~ class, data = survival)
      pdf(paste(cancer, "_senescence_sur.pdf", sep = ""), height = 6, width = 8, onefile = FALSE)
      p5 <- ggsurvplot(
        fit,
        data = survival,
        palette = c("#E7B800", "#2E9FDF"),
        conf.int = FALSE,
        pval = T,
        pval.method = T,
        risk.table = T,
        cumevents = F,
        cumcensor = F,
      )
      print(p5)
      dev.off()
    }
  }
}

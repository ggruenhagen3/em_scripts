#*******************************************************************************
# Load Libraries ===============================================================
#*******************************************************************************
wdstr = substr(getwd(), 1, 12)
switch(wdstr,
       "C:/Users/mil" = { main_path = "C:/Users/miles/Downloads/";        },
       "/home/george" = { main_path = "~/research/"                       },
       "/storage/scr" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/hom" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/cod" = { main_path = "/storage/scratch1/6/ggruenhagen3/" })
brain_dir = paste0(main_path, "brain/")
data_dir  = paste0(main_path, "em/data/")
out_dir   = paste0(main_path, "em/results/")
if (main_path == "/storage/scratch1/6/ggruenhagen3/") { data_dir = "/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/data/st/data/" }
source(paste0(brain_dir, "/brain_scripts/all_f.R"))
setwd(out_dir)

#*******************************************************************************
# Load Objects =================================================================
#*******************************************************************************
gene_info = read.table(paste0(main_path, "/all_research/gene_info_2.txt"), header = T, stringsAsFactors = F)
em = readRDS(paste0(data_dir, "em_120822.rds"))

#*******************************************************************************
# Initial Clustering ===========================================================
#*******************************************************************************
dir_of_sr_dirs = "~/scratch/brain/bs/JTS14/" # Folder where all the individual samples are kept
em.counts = Read10X(paste0(dir_of_sr_dirs, "/PM19_B1_nuc/outs/filtered_feature_bc_matrix/"))
em = CreateSeuratObject(em.counts)

#Check the quality
plot1 = ggplot(data.frame(nCount = em$nCount_RNA), aes(x = nCount)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nCount")
plot2 = VlnPlot(em, features = "nCount_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())

plot3 = ggplot(data.frame(nFeature = em$nFeature_RNA), aes(x = nFeature)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nFeature")
plot4 = VlnPlot(em, features = "nFeature_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
pdf(paste0("~/scratch/brain/results/hb_", i, "quality.pdf"), width = 8, height = 8)
print(cowplot::plot_grid(plotlist = list(plot1, plot2, plot3, plot4), ncol = 2))
dev.off()
print(paste0("# Cells UMI > 300 = ", length(which(em$nCount_RNA > 300))))
print(paste0("# Cells Genes > 300 = ", length(which(em$nFeature_RNA > 300))))

gtf = read.delim("~/scratch/m_zebra_ref/gtf_processed.gtf", header = F)
mito.genes = gtf$V10[which(gtf$V1 == "NC_027944.1")]
mito.genes = stringr::str_replace(mito.genes, "_", "-")
em$pct.mt = colSums(em@assays$RNA@counts[mito.genes,]) / em$nCount_RNA

# Clustering
em = subset(em, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & pct.mt < 0.05)
print(paste("Number of Cells in em After Filterning:", ncol(em)))
em = NormalizeData(em, normalization.method = "LogNormalize", scale.factor = 10000)
em = SCTransform(em, vars.to.regress = "sample" , verbose = TRUE)
em@active.assay = "SCT"
em = RunPCA(em)
em = RunUMAP(em, reduction = "pca", dims = 1:30)
em = FindNeigemors(em, reduction="umap", dims = 1:2)
em = FindClusters(em, resolution = 0.30, algorithm = 2)

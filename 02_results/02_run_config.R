source("/home/claude/02_cna_eligibility_scan.R")
base <- "/mnt/user-data/uploads/KIDS26-Team13/examples"
cfg <- data.frame(
  gene = c("PTEN","CDKN2A"),
  transcript = c("ENST00000371953","ENST00000304494"),
  exon.file  = file.path(base, c("PTEN_CNA/input_files/PTEN_exon_to_protein_map.csv",
                                 "CDKN2A_CNA/input_files/CDKN2A_exon_to_protein_map.csv")),
  coord.file = file.path(base, c("PTEN_CNA/input_files/PTEN_alphafold_coordinates.csv",
                                 "CDKN2A_CNA/input_files/CDKN2A_alphafold_coordinates.csv")),
  cna.file   = file.path(base, c("PTEN_CNA/input_files/PTEN_CNA_lesions.csv",
                                 "CDKN2A_CNA/input_files/CDKN2A_CNA_lesions.csv")),
  stringsAsFactors = FALSE)
cols <- list(subject="ID", event=NULL, chrom="chrom", start="loc.start", end="loc.end", type="lsn.type")
# grid extended: min.plddt = 0 reproduces the shipped examples (use.plddt.filter = FALSE)
grid <- expand.grid(min.plddt=c(0,50,70,90),
                    min.exon.coding.overlap=c(0,0.10,0.50),
                    min.eligible.exon.fraction=c(1.00,0.80),
                    min.exon.mapped.fraction=c(0.50), stringsAsFactors=FALSE)
p <- baseline.params; p$min.plddt <- 0   # baseline = the example's actual setting
res <- run_eligibility_scan(genes=NULL, cna.columns=cols, gene.config=cfg,
                            params=p, grid=grid, results.dir="/home/claude/02_results")
saveRDS(res, "/home/claude/res02.rds")
print(t(res$overview))

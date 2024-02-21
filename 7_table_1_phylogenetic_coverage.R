source("~/sync/mangroves/M2_Alice/resultats/donnees_scripts/github/nanopore_metabarcoding/5_main_phyloseq_objects.R")

#########################################################################@
########################################@
## Phylogenetic coverage of bacteria from Illumina and Nanopore sequencing
## Nanopore detected 95.6% of the phyla detected by Illumina
## Nanopore detected 93.4% of the classes detected by Illumina
## Nanopore detected 92.8% of the orders detected by Illumina
## Nanopore detected 92.2% of the families detected by Illumina
## Nanopore detected 87.7% of the genus detected by Illumina
## Nanopore detected 69.7% of the species detected by Illumina



#### SPECIES (bacteria only)
species_illu <- paste(illu_bact_raref_2reads@tax_table[,2],illu_bact_raref_2reads@tax_table[,3],illu_bact_raref_2reads@tax_table[,4],
                      illu_bact_raref_2reads@tax_table[,5],illu_bact_raref_2reads@tax_table[,6],illu_bact_raref_2reads@tax_table[,7])
species_nano <- paste(nano_bact_raref_2reads@tax_table[,2],nano_bact_raref_2reads@tax_table[,3],nano_bact_raref_2reads@tax_table[,4],
                      nano_bact_raref_2reads@tax_table[,5],nano_bact_raref_2reads@tax_table[,6],nano_bact_raref_2reads@tax_table[,7])
species_all <- paste(bact_raref_2reads@tax_table[,2],bact_raref_2reads@tax_table[,3],bact_raref_2reads@tax_table[,4],
                     bact_raref_2reads@tax_table[,5],bact_raref_2reads@tax_table[,6],bact_raref_2reads@tax_table[,7])

n_occur <- data.frame(table(species_all))
n_duplicates <- n_occur[n_occur$Freq > 1,] ## no duplicates

taxa_names(illu_bact_raref_2reads) <- species_illu
taxa_names(nano_bact_raref_2reads) <- species_nano
taxa_names(bact_raref_2reads) <- species_all

species_illu %>% length() ## 749
species_nano %>% length() ## 1495
species_all %>% length() ## 2244
setdiff(species_illu,species_nano) %>% length() # 227 species detected by illu only
setdiff(species_nano,species_illu) %>% length() # 973 species detected by nano only
intersect(species_illu,species_nano) %>% length() ## 522

(749-227)/749 ## Nanopore detected 69.7% of the species detected by Illumina
749-522
1495-522

sum(taxa_sums(subset_taxa(nano_bact_raref_2reads,Species!="__") %>% tax_select(tax_list = intersect(species_illu,species_nano), n_typos = 0, ranks_searched = "Species"))) /
  sum(taxa_sums(subset_taxa(nano_bact_raref_2reads,Species!="__"))) # 68.5% of the Nanopore reads are assigned to species detected by Illumina (among assigned species)
sum(taxa_sums(subset_taxa(illu_bact_raref_2reads,Species!="__") %>% tax_select(tax_list = intersect(species_illu,species_nano), n_typos = 0, ranks_searched = "Species"))) /
  sum(taxa_sums(subset_taxa(illu_bact_raref_2reads,Species!="__"))) # 84.7% of the Illumina reads are assigned to species detected by Nanopore (among assigned species)

sum(taxa_sums(subset_taxa(nano_bact_raref_2reads,Phylum!="__") %>% tax_select(tax_list = intersect(species_illu,species_nano), n_typos = 0, ranks_searched = "Phylum"))) /
  sum(taxa_sums(subset_taxa(nano_bact_raref_2reads,Phylum!="__"))) # 75.2% of the Nanopore reads are assigned to Phylum detected by Illumina (among assigned Phyla)
sum(taxa_sums(subset_taxa(illu_bact_raref_2reads,Phylum!="__") %>% tax_select(tax_list = intersect(species_illu,species_nano), n_typos = 0, ranks_searched = "Phylum"))) /
  sum(taxa_sums(subset_taxa(illu_bact_raref_2reads,Phylum!="__"))) # 90.3% of the Illumina reads are assigned to Phylum detected by Nanopore (among assigned Phyla)


library(stringr)
table(str_count(setdiff(species_nano,species_illu), "__"))
## among the 973 species detected by nano only, 892 (91.7%) were identified to the genus rank.
660+232+47+16+11+7
660/973
892/973
7/11
11/41
660/973
table(str_count(setdiff(species_illu,species_nano), "__"))
table(str_count(intersect(species_illu,species_nano), "__"))
1-306/522
1-147/393
1-35/285

# For nanopore, 529 unknown species over 1495, i.e. 35.4% unassigned species
1 - ((subset_taxa(nano_bact_raref_2reads,Species=="__") %>% taxa_names() %>% length())/ 
       (nano_bact_raref_2reads %>% taxa_names() %>% length())) # 64.6% species assigned
1 - (subset_taxa(nano_bact_raref_2reads,Species=="__") %>% sample_sums() %>% sum())/
  (nano_bact_raref_2reads %>% sample_sums() %>% sum()) # 46.9% reads assigned
# For illumina, 260 unknown species over 749, i.e. 34.7% unassigned species
1 - ((subset_taxa(illu_bact_raref_2reads,Species=="__") %>% taxa_names() %>% length())/ 
       (illu_bact_raref_2reads %>% taxa_names() %>% length())) # 65.3% species assigned
1 - (subset_taxa(illu_bact_raref_2reads,Species=="__") %>% sample_sums() %>% sum())/
  (illu_bact_raref_2reads %>% sample_sums() %>% sum()) # 54.0% reads assigned

## shared / unshared species
illu_bact_raref_2reads_sh <- tax_select(illu_bact_raref_2reads, intersect(species_illu,species_nano))
nano_bact_raref_2reads_sh <- tax_select(nano_bact_raref_2reads, intersect(species_nano,species_illu))
illu_bact_raref_2reads_unsh <- tax_select(illu_bact_raref_2reads, setdiff(species_illu,species_nano))
nano_bact_raref_2reads_unsh <- tax_select(nano_bact_raref_2reads, setdiff(species_nano,species_illu))
1 - ((subset_taxa(illu_bact_raref_2reads_sh,Species=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_sh %>% taxa_names() %>% length())) # Among the 522 shared species, 306 (58.6%) were assigned by Illumina, 
306/522
1 - ((subset_taxa(illu_bact_raref_2reads_unsh,Species=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_unsh %>% taxa_names() %>% length())) #
1 - ((subset_taxa(illu_bact_raref_2reads_sh,Species=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_sh %>% sample_sums() %>% sum())) # 50.6% of Illumina reads assigned for shared taxa
1 - ((subset_taxa(nano_bact_raref_2reads_sh,Species=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_sh %>% sample_sums() %>% sum())) # 41.1% of Nanopore reads assigned for shared taxa
1 - ((subset_taxa(nano_bact_raref_2reads_unsh,Species=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_unsh %>% taxa_names() %>% length())) #
1 - ((subset_taxa(illu_bact_raref_2reads_unsh,Species=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_unsh %>% sample_sums() %>% sum())) # 85.8% of Illumina reads assigned for unshared taxa
1 - ((subset_taxa(nano_bact_raref_2reads_unsh,Species=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_unsh %>% sample_sums() %>% sum())) # 67.4% of Nanopore reads assigned for unshared taxa


#### GENUS
illu_bact_raref_2reads_gen <- tax_glom(illu_bact_raref_2reads, taxrank = "Genus")
nano_bact_raref_2reads_gen <- tax_glom(nano_bact_raref_2reads, taxrank = "Genus")

genus_illu <- paste(illu_bact_raref_2reads_gen@tax_table[,2],illu_bact_raref_2reads_gen@tax_table[,3],illu_bact_raref_2reads_gen@tax_table[,4],
                    illu_bact_raref_2reads_gen@tax_table[,5],illu_bact_raref_2reads_gen@tax_table[,6])
genus_nano <- paste(nano_bact_raref_2reads_gen@tax_table[,2],nano_bact_raref_2reads_gen@tax_table[,3],nano_bact_raref_2reads_gen@tax_table[,4],
                    nano_bact_raref_2reads_gen@tax_table[,5],nano_bact_raref_2reads_gen@tax_table[,6])
taxa_names(illu_bact_raref_2reads_gen) <- genus_illu
taxa_names(nano_bact_raref_2reads_gen) <- genus_nano
genus_illu %>% length() ## 448
genus_nano %>% length() ## 785
setdiff(genus_illu,genus_nano) %>% length() # 55 genus detected by illu only
setdiff(genus_nano,genus_illu) %>% length() # 392 genus detected by nano only
(448-55)/448
(448-392) + 785

## Nanopore detected 87.7% of the genus detected by Illumina

sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_gen,Genus!="__") %>% tax_select(tax_list = intersect(genus_illu,genus_nano), n_typos = 0, ranks_searched = "Genus"))) /
  sum(taxa_sums(subset_taxa(nano_bact_raref_2reads,Genus!="__"))) # 89.6% of the Nanopore reads are assigned to Genus detected by Illumina (among assigned genus)
sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_gen,Genus!="__") %>% tax_select(tax_list = intersect(genus_illu,genus_nano), n_typos = 0, ranks_searched = "Genus"))) /
  sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_gen,Genus!="__"))) # 98.8% of the Illumina reads are assigned to Genus detected by Nanopore (among assigned genus)


1 - ((subset_taxa(nano_bact_raref_2reads_gen,Genus=="__") %>% taxa_names() %>% length())/ 
       (nano_bact_raref_2reads_gen %>% taxa_names() %>% length())) # 80.9% genus assigned
1 - (subset_taxa(nano_bact_raref_2reads_gen,Genus=="__") %>% sample_sums() %>% sum())/
  (nano_bact_raref_2reads_gen %>% sample_sums() %>% sum()) # 67.8% reads assigned
1 - ((subset_taxa(illu_bact_raref_2reads_gen,Genus=="__") %>% taxa_names() %>% length())/ 
       (illu_bact_raref_2reads_gen %>% taxa_names() %>% length())) # 82.4% genus assigned
1 - (subset_taxa(illu_bact_raref_2reads_gen,Genus=="__") %>% sample_sums() %>% sum())/
  (illu_bact_raref_2reads_gen %>% sample_sums() %>% sum()) # 94.0% reads assigned

## shared / unshared genus
illu_bact_raref_2reads_gen_sh <- tax_select(illu_bact_raref_2reads_gen, intersect(genus_illu,genus_nano))
nano_bact_raref_2reads_gen_sh <- tax_select(nano_bact_raref_2reads_gen, intersect(genus_nano,genus_illu))
illu_bact_raref_2reads_gen_unsh <- tax_select(illu_bact_raref_2reads_gen, setdiff(genus_illu,genus_nano))
nano_bact_raref_2reads_gen_unsh <- tax_select(nano_bact_raref_2reads_gen, setdiff(genus_nano,genus_illu))
1 - ((subset_taxa(illu_bact_raref_2reads_gen_sh,Genus=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_gen_sh %>% taxa_names() %>% length())) # Among the 393 shared genus, 379 (82.4%) were assigned by Illumina, 
1 - ((subset_taxa(nano_bact_raref_2reads_gen_sh,Genus=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_gen_sh %>% taxa_names() %>% length())) # idem for Nanopore shared genus
illu_bact_raref_2reads_gen_unsh %>% taxa_names() %>% length() # 55 genus unshared for Illumina
nano_bact_raref_2reads_gen_unsh %>% taxa_names() %>% length() # 392 genus unshared for Nanopore
(448-69) # ?
1 - ((subset_taxa(illu_bact_raref_2reads_gen_unsh,Genus=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_gen_unsh %>% taxa_names() %>% length())) # 81.8% of Illumina taxa assigned for unshared genus
1 - ((subset_taxa(illu_bact_raref_2reads_gen_sh,Genus=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_gen_sh %>% sample_sums() %>% sum())) # 94.2% of Illumina reads assigned for shared genus
1 - ((subset_taxa(nano_bact_raref_2reads_gen_unsh,Genus=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_gen_unsh %>% taxa_names() %>% length())) # 79.3% of Nanopore taxa assigned for unshared genus
1 - ((subset_taxa(nano_bact_raref_2reads_gen_sh,Genus=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_gen_sh %>% sample_sums() %>% sum())) # 66.5% of Nanopore reads assigned for shared genus
1 - ((subset_taxa(illu_bact_raref_2reads_gen_unsh,Genus=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_gen_unsh %>% sample_sums() %>% sum())) # 84.5% of Illumina reads assigned for unshared genus
1 - ((subset_taxa(nano_bact_raref_2reads_gen_unsh,Genus=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_gen_unsh %>% sample_sums() %>% sum())) # 81.8% of Nanopore reads assigned for unshared genus



#### FAMILIES
illu_bact_raref_2reads_fam <- tax_glom(illu_bact_raref_2reads, taxrank = "Family")
nano_bact_raref_2reads_fam <- tax_glom(nano_bact_raref_2reads, taxrank = "Family")
fam_illu <- paste(illu_bact_raref_2reads_fam@tax_table[,2],illu_bact_raref_2reads_fam@tax_table[,3],illu_bact_raref_2reads_fam@tax_table[,4],
                  illu_bact_raref_2reads_fam@tax_table[,5])
fam_nano <- paste(nano_bact_raref_2reads_fam@tax_table[,2],nano_bact_raref_2reads_fam@tax_table[,3],nano_bact_raref_2reads_fam@tax_table[,4],
                  nano_bact_raref_2reads_fam@tax_table[,5])
taxa_names(illu_bact_raref_2reads_fam) <- fam_illu
taxa_names(nano_bact_raref_2reads_fam) <- fam_nano
fam_illu %>% length() ## 309
fam_nano %>% length() ## 483
setdiff(fam_illu,fam_nano) # 24 families detected by illu only
setdiff(fam_nano,fam_illu) # 198 families detected by nano only
(309-24)/309
## Nanopore detected 92.2% of the families detected by Illumina

sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_fam,Family!="__") %>% tax_select(tax_list = intersect(fam_illu,fam_nano), n_typos = 0, ranks_searched = "Family"))) /
  sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_fam,Family!="__"))) # 94.7 of the Nanopore reads are assigned to Family detected by Illumina (among assigned Family)
sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_fam,Family!="__") %>% tax_select(tax_list = intersect(fam_illu,fam_nano), n_typos = 0, ranks_searched = "Family"))) /
  sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_fam,Family!="__"))) # 99.7% of the Illumina reads are assigned to Family detected by Nanopore (among assigned Family)


1 - ((subset_taxa(nano_bact_raref_2reads_fam,Family=="__") %>% taxa_names() %>% length())/ 
       (nano_bact_raref_2reads_fam %>% taxa_names() %>% length())) # 85.9% families assigned
1 - (subset_taxa(nano_bact_raref_2reads_fam,Family=="__") %>% sample_sums() %>% sum())/
  (nano_bact_raref_2reads_fam %>% sample_sums() %>% sum()) # 73.3% reads assigned
1 - ((subset_taxa(illu_bact_raref_2reads_fam,Family=="__") %>% taxa_names() %>% length())/ 
       (illu_bact_raref_2reads_fam %>% taxa_names() %>% length())) # 88.0% familes assigned
1 - (subset_taxa(illu_bact_raref_2reads_fam,Family=="__") %>% sample_sums() %>% sum())/
  (illu_bact_raref_2reads_fam %>% sample_sums() %>% sum()) # 98.0% reads assigned

## shared / unshared families
illu_bact_raref_2reads_fam_sh <- tax_select(illu_bact_raref_2reads_fam, intersect(fam_illu,fam_nano))
nano_bact_raref_2reads_fam_sh <- tax_select(nano_bact_raref_2reads_fam, intersect(fam_nano,fam_illu))
illu_bact_raref_2reads_fam_unsh <- tax_select(illu_bact_raref_2reads_fam, setdiff(fam_illu,fam_nano))
nano_bact_raref_2reads_fam_unsh <- tax_select(nano_bact_raref_2reads_fam, setdiff(fam_nano,fam_illu))
1 - ((subset_taxa(illu_bact_raref_2reads_fam_sh,Family=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_fam_sh %>% taxa_names() %>% length())) # Among the 285 shared Family, 251 (88.1%) were assigned by Illumina, 
1 - ((subset_taxa(nano_bact_raref_2reads_fam_sh,Family=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_fam_sh %>% taxa_names() %>% length())) # idem for Nanopore shared Family
illu_bact_raref_2reads_fam_unsh %>% taxa_names() %>% length() # 24 Family unshared for Illumina
nano_bact_raref_2reads_fam_unsh %>% taxa_names() %>% length() # 198 Family unshared for Nanopore
(285-34) # ?
1 - ((subset_taxa(illu_bact_raref_2reads_fam_unsh,Family=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_fam_unsh %>% taxa_names() %>% length())) # 87.5% of Illumina taxa assigned for unshared Family
1 - ((subset_taxa(illu_bact_raref_2reads_fam_sh,Family=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_fam_sh %>% sample_sums() %>% sum())) # 98.0% of Illumina reads assigned for shared Family
1 - ((subset_taxa(nano_bact_raref_2reads_fam_unsh,Family=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_fam_unsh %>% taxa_names() %>% length())) # 82.8% of Nanopore taxa assigned for unshared Family
1 - ((subset_taxa(nano_bact_raref_2reads_fam_sh,Family=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_fam_sh %>% sample_sums() %>% sum())) # 72.8% of Nanopore reads assigned for shared Family
1 - ((subset_taxa(illu_bact_raref_2reads_fam_unsh,Family=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_fam_unsh %>% sample_sums() %>% sum())) # 83.9% of Illumina reads assigned for unshared Family
1 - ((subset_taxa(nano_bact_raref_2reads_fam_unsh,Family=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_fam_unsh %>% sample_sums() %>% sum())) # 83.7% of Nanopore reads assigned for unshared Family



#### ORDERS
illu_bact_raref_2reads_ord <- tax_glom(illu_bact_raref_2reads, taxrank = "Order")
nano_bact_raref_2reads_ord <- tax_glom(nano_bact_raref_2reads, taxrank = "Order")
orders_illu <- paste(illu_bact_raref_2reads_ord@tax_table[,2],illu_bact_raref_2reads_ord@tax_table[,3],illu_bact_raref_2reads_ord@tax_table[,4])
orders_nano <- paste(nano_bact_raref_2reads_ord@tax_table[,2],nano_bact_raref_2reads_ord@tax_table[,3],nano_bact_raref_2reads_ord@tax_table[,4])
orders_illu %>% length() ## 209
orders_nano %>% length() ## 316
taxa_names(illu_bact_raref_2reads_ord) <- orders_illu
taxa_names(nano_bact_raref_2reads_ord) <- orders_nano
setdiff(orders_illu,orders_nano) %>% length() # 15 orders detected by illu only
setdiff(orders_nano,orders_illu) %>% length() # 122 orders detected by nano only
(209-15)/209
## Nanopore detected 92.8% of the orders detected by Illumina

sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_ord,Order!="__") %>% tax_select(tax_list = intersect(orders_illu,orders_nano), n_typos = 0, ranks_searched = "Order"))) /
  sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_ord,Order!="__"))) # 95.7 of the Nanopore reads are assigned to Order detected by Illumina (among assigned Order)
sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_ord,Order!="__") %>% tax_select(tax_list = intersect(orders_illu,orders_nano), n_typos = 0, ranks_searched = "Order"))) /
  sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_ord,Order!="__"))) # 99.7% of the Illumina reads are assigned to Order detected by Nanopore (among assigned Order)



1 - ((subset_taxa(nano_bact_raref_2reads_ord,Order=="__") %>% taxa_names() %>% length())/ 
       (nano_bact_raref_2reads_ord %>% taxa_names() %>% length())) # 89.2% orders assigned
1 - (subset_taxa(nano_bact_raref_2reads_ord,Order=="__") %>% sample_sums() %>% sum())/
  (nano_bact_raref_2reads_ord %>% sample_sums() %>% sum()) # 74.9% reads assigned
1 - ((subset_taxa(illu_bact_raref_2reads_ord,Order=="__") %>% taxa_names() %>% length())/ 
       (illu_bact_raref_2reads_ord %>% taxa_names() %>% length())) # 91.9% orders assigned
1 - (subset_taxa(illu_bact_raref_2reads_ord,Order=="__") %>% sample_sums() %>% sum())/
  (illu_bact_raref_2reads_ord %>% sample_sums() %>% sum()) # 98.3% reads assigned

## shared / unshared orders
illu_bact_raref_2reads_ord_sh <- tax_select(illu_bact_raref_2reads_ord, intersect(orders_illu,orders_nano), ranks_searched="Class", n_typos = 0)
nano_bact_raref_2reads_ord_sh <- tax_select(nano_bact_raref_2reads_ord, intersect(orders_nano,orders_illu), ranks_searched="Class", n_typos = 0)
illu_bact_raref_2reads_ord_unsh <- tax_select(illu_bact_raref_2reads_ord, setdiff(orders_illu,orders_nano), ranks_searched="Class", n_typos = 0)
nano_bact_raref_2reads_ord_unsh <- tax_select(nano_bact_raref_2reads_ord, setdiff(orders_nano,orders_illu), ranks_searched="Class", n_typos = 0)
1 - ((subset_taxa(illu_bact_raref_2reads_ord_sh,Order=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_ord_sh %>% taxa_names() %>% length())) # Among the 194 shared Orders, 178 (91.8%) were assigned by Illumina, 
1 - ((subset_taxa(nano_bact_raref_2reads_ord_sh,Order=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_ord_sh %>% taxa_names() %>% length())) # idem for Nanopore shared Order
illu_bact_raref_2reads_ord_unsh %>% taxa_names() %>% length() # 15 Orders unshared for Illumina
nano_bact_raref_2reads_ord_unsh %>% taxa_names() %>% length() # 122 Orders unshared for Nanopore
(194-16) # ?
1 - ((subset_taxa(illu_bact_raref_2reads_ord_unsh,Order=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_ord_unsh %>% taxa_names() %>% length())) # 93.3% of Illumina taxa assigned for unshared Orders
1 - ((subset_taxa(illu_bact_raref_2reads_ord_sh,Order=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_ord_sh %>% sample_sums() %>% sum())) # 98.3% of Illumina reads assigned for shared Orders
1 - ((subset_taxa(nano_bact_raref_2reads_ord_unsh,Order=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_ord_unsh %>% taxa_names() %>% length())) # 85.2% of Nanopore taxa assigned for unshared Orders
1 - ((subset_taxa(nano_bact_raref_2reads_ord_sh,Order=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_ord_sh %>% sample_sums() %>% sum())) # 74.4% of Nanopore reads assigned for shared Orders
1 - ((subset_taxa(illu_bact_raref_2reads_ord_unsh,Order=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_ord_unsh %>% sample_sums() %>% sum())) # 97.1% of Illumina reads assigned for unshared Orders
1 - ((subset_taxa(nano_bact_raref_2reads_ord_unsh,Order=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_ord_unsh %>% sample_sums() %>% sum())) # 87.2% of Nanopore reads assigned for unshared Orders



#### CLASSES
illu_bact_raref_2reads_class <- tax_glom(illu_bact_raref_2reads, taxrank = "Class", NArm=FALSE)
nano_bact_raref_2reads_class <- tax_glom(nano_bact_raref_2reads, taxrank = "Class", NArm=FALSE)
class_illu <- paste(illu_bact_raref_2reads_class@tax_table[,2],illu_bact_raref_2reads_class@tax_table[,3])
class_nano <- paste(nano_bact_raref_2reads_class@tax_table[,2],nano_bact_raref_2reads_class@tax_table[,3])
class_illu %>% length() ## 106
class_nano %>% length() ## 140
taxa_names(illu_bact_raref_2reads_class) <- class_illu
taxa_names(nano_bact_raref_2reads_class) <- class_nano
setdiff(class_illu,class_nano) %>% length() # 7 classes detected by illu only
setdiff(class_nano,class_illu) %>% length() # 41 classes detected by nano only
(106-7)/106
## Nanopore detected 93.4% of the classes detected by Illumina

sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_class,Class!="__") %>% tax_select(tax_list = intersect(class_illu,class_nano), n_typos = 0, ranks_searched = "Class"))) /
  sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_class,Class!="__"))) # 99.3% of the Nanopore reads are assigned to Class detected by Illumina (among assigned Class)
sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_class,Class!="__") %>% tax_select(tax_list = intersect(class_illu,class_nano), n_typos = 0, ranks_searched = "Class"))) /
  sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_class,Class!="__"))) # 99.97% of the Illumina reads are assigned to Class detected by Nanopore (among assigned Class)


1 - ((subset_taxa(nano_bact_raref_2reads_class,Class=="__") %>% taxa_names() %>% length())/ 
       (nano_bact_raref_2reads_class %>% taxa_names() %>% length())) # 89.3% classes assigned
1 - (subset_taxa(nano_bact_raref_2reads_class,Class=="__") %>% sample_sums() %>% sum())/
  (nano_bact_raref_2reads_class %>% sample_sums() %>% sum()) # 86.0% reads assigned
1 - ((subset_taxa(illu_bact_raref_2reads_class,Class=="__") %>% taxa_names() %>% length())/ 
       (illu_bact_raref_2reads_class %>% taxa_names() %>% length())) # 92.5% classes assigned
1 - (subset_taxa(illu_bact_raref_2reads_class,Class=="__") %>% sample_sums() %>% sum())/
  (illu_bact_raref_2reads_class %>% sample_sums() %>% sum()) # 99.6% reads assigned

## shared / unshared classes
illu_bact_raref_2reads_class_sh <- tax_select(illu_bact_raref_2reads_class, intersect(class_illu,class_nano), ranks_searched="Class", n_typos = 0, strict_matches=TRUE)
nano_bact_raref_2reads_class_sh <- tax_select(nano_bact_raref_2reads_class, intersect(class_illu,class_nano), ranks_searched="Class", n_typos = 0, strict_matches=TRUE)
illu_bact_raref_2reads_class_unsh <- tax_select(illu_bact_raref_2reads_class, setdiff(class_illu,class_nano), ranks_searched="Class", n_typos = 0, strict_matches=TRUE)
nano_bact_raref_2reads_class_unsh <- tax_select(nano_bact_raref_2reads_class, setdiff(class_nano,class_illu), ranks_searched="Class", n_typos = 0, strict_matches=TRUE)
1 - ((subset_taxa(illu_bact_raref_2reads_class_sh,Class=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_class_sh %>% taxa_names() %>% length())) # Among the 99 shared class, 91 (91.9%) were assigned by Illumina, 
1 - ((subset_taxa(nano_bact_raref_2reads_class_sh,Class=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_class_sh %>% taxa_names() %>% length())) # idem for Nanopore shared Class
illu_bact_raref_2reads_class_unsh %>% taxa_names() %>% length() # 7 class unshared for Illumina
nano_bact_raref_2reads_class_unsh %>% taxa_names() %>% length() # 41 class unshared for Nanopore
1 - ((subset_taxa(illu_bact_raref_2reads_class_unsh,Class=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_class_unsh %>% taxa_names() %>% length())) # 100% of Illumina taxa assigned for unshared class
1 - ((subset_taxa(illu_bact_raref_2reads_class_unsh,Class=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_class_unsh %>% sample_sums() %>% sum())) # 100% of Illumina reads assigned for unshared class
1 - ((subset_taxa(illu_bact_raref_2reads_class_sh,Class=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_class_sh %>% sample_sums() %>% sum())) # 99.6% of Illumina reads assigned for shared class
1 - ((subset_taxa(nano_bact_raref_2reads_class_unsh,Class=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_class_unsh %>% taxa_names() %>% length())) # 82.9% of Nanopore taxa assigned for unshared class
1 - ((subset_taxa(nano_bact_raref_2reads_class_unsh,Class=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_class_unsh %>% sample_sums() %>% sum())) # 70.5% of Nanopore reads assigned for unshared class
1 - ((subset_taxa(nano_bact_raref_2reads_class_sh,Class=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_class_sh %>% sample_sums() %>% sum())) # 86.1% of Nanopore reads assigned for shared class

nano_bact_raref_2reads_class
nano_bact_raref_2reads_class %>% sample_sums() %>% sum()
((nano_bact_raref_2reads_class_unsh %>% sample_sums() %>% sum()) + (nano_bact_raref_2reads_class_sh %>% sample_sums() %>% sum()))/(nano_bact_raref_2reads_class %>% sample_sums() %>% sum())
(nano_bact_raref_2reads_ord_unsh %>% sample_sums() %>% sum()) + (nano_bact_raref_2reads_ord_sh %>% sample_sums() %>% sum())

subset_taxa(nano_bact_raref_2reads_class_unsh,Class=="__") %>% taxa_names()
subset_taxa(nano_bact_raref_2reads_ord_unsh, Order=="__") %>% taxa_names() %>% substr(., 1, nchar(.)-3)


#### PHYLA
illu_bact_raref_2reads_phy <- tax_glom(illu_bact_raref_2reads, taxrank = "Phylum")
nano_bact_raref_2reads_phy <- tax_glom(nano_bact_raref_2reads, taxrank = "Phylum")
phyla_illu <- illu_bact_raref_2reads_phy@tax_table[,2]
phyla_nano <- nano_bact_raref_2reads_phy@tax_table[,2]
phyla_illu %>% length() ## 45
phyla_nano %>% length() ## 54
taxa_names(illu_bact_raref_2reads_phy) <- phyla_illu
taxa_names(nano_bact_raref_2reads_phy) <- phyla_nano
setdiff(phyla_illu,phyla_nano) %>% length() # 2 phyla detected by illu only (Cloacimonadota, CK-2C2-2)
setdiff(phyla_nano,phyla_illu) %>% length() # 11 phyla detected by nano only Acetothermia, WS2, LCP-89, WOR-1, Armatimonadota, Margulisbacteria,
## Nitrospinota, Fermentibacterota, Methylomirabilota, Caldatribacteriota, WPS-2
(45-2)/45
## Nanopore detected 95.6% of the phyla detected by Illumina

sum(taxa_sums(illu_bact_raref_2reads_phy))
sum(taxa_sums(nano_bact_raref_2reads_phy))
sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_phy,Phylum!="__") %>% tax_select(tax_list = intersect(phyla_illu,phyla_nano), n_typos = 0, ranks_searched = "Phylum"))) /
  sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_phy,Phylum!="__"))) # 99.7% of the Nanopore reads are assigned to Phylum detected by Illumina (among assigned Phylum)
sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_phy,Phylum!="__") %>% tax_select(tax_list = intersect(phyla_illu,phyla_nano), n_typos = 0, ranks_searched = "Phylum"))) /
  sum(taxa_sums(subset_taxa(illu_bact_raref_2reads_phy,Phylum!="__"))) # 99.99% of the Illumina reads are assigned to Phylum detected by Nanopore (among assigned Phylum)



1 - ((subset_taxa(nano_bact_raref_2reads_phy,Phylum=="__") %>% taxa_names() %>% length())/ 
       (nano_bact_raref_2reads_phy %>% taxa_names() %>% length())) # 98.1% Phylums assigned
1 - (subset_taxa(nano_bact_raref_2reads_phy,Phylum=="__") %>% sample_sums() %>% sum())/
  (nano_bact_raref_2reads_phy %>% sample_sums() %>% sum()) # 88.3% reads assigned
1 - ((subset_taxa(illu_bact_raref_2reads_phy,Phylum=="__") %>% taxa_names() %>% length())/ 
       (illu_bact_raref_2reads_phy %>% taxa_names() %>% length())) # 97.8% Phylums assigned
1 - (subset_taxa(illu_bact_raref_2reads_phy,Phylum=="__") %>% sample_sums() %>% sum())/
  (illu_bact_raref_2reads_phy %>% sample_sums() %>% sum()) # 99.6% reads assigned

## shared / unshared phyla
illu_bact_raref_2reads_phy_sh <- tax_select(illu_bact_raref_2reads_phy, intersect(phyla_illu,phyla_nano), ranks_searched="Phylum", n_typos = 0)
nano_bact_raref_2reads_phy_sh <- tax_select(nano_bact_raref_2reads_phy, intersect(phyla_illu,phyla_nano), ranks_searched="Phylum", n_typos = 0)
illu_bact_raref_2reads_phy_unsh <- tax_select(illu_bact_raref_2reads_phy, setdiff(phyla_illu,phyla_nano), ranks_searched="Phylum", n_typos = 0)
nano_bact_raref_2reads_phy_unsh <- tax_select(nano_bact_raref_2reads_phy, setdiff(phyla_nano,phyla_illu), ranks_searched="Phylum", n_typos = 0) # 13 phyla au lieu de 11 ?

intersect(phyla_nano,phyla_illu) %>% length()
setdiff((nano_bact_raref_2reads_phy_sh %>% taxa_names()), (illu_bact_raref_2reads_phy_sh %>% taxa_names()))

1 - ((subset_taxa(illu_bact_raref_2reads_phy_sh,Phylum=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_phy_sh %>% taxa_names() %>% length())) # Among the 43 shared phyla, 97.6% were assigned by Illumina, 
1 - ((subset_taxa(nano_bact_raref_2reads_phy_sh,Phylum=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_phy_sh %>% taxa_names() %>% length())) # idem for Nanopore shared Phyla
illu_bact_raref_2reads_phy_unsh %>% taxa_names() %>% length() # 2 phyla unshared for Illumina
nano_bact_raref_2reads_phy_unsh %>% taxa_names() %>% length() # 11 phyla unshared for Nanopore
1 - ((subset_taxa(illu_bact_raref_2reads_phy_unsh,Phylum=="__") %>% taxa_names() %>% length())/
       (illu_bact_raref_2reads_phy_unsh %>% taxa_names() %>% length())) # 100% of Illumina taxa assigned for unshared Phyla
1 - ((subset_taxa(illu_bact_raref_2reads_phy_unsh,Phylum=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_phy_unsh %>% sample_sums() %>% sum())) # 100% of Illumina reads assigned for unshared Phyla
1 - ((subset_taxa(illu_bact_raref_2reads_phy_sh,Phylum=="__")) %>% sample_sums() %>% sum() /
       (illu_bact_raref_2reads_phy_sh %>% sample_sums() %>% sum())) # 99.6% of Illumina reads assigned for shared Phyla
1 - ((subset_taxa(nano_bact_raref_2reads_phy_unsh,Phylum=="__") %>% taxa_names() %>% length())/
       (nano_bact_raref_2reads_phy_unsh %>% taxa_names() %>% length())) # 100% of Nanopore taxa assigned for unshared Phyla
1 - ((subset_taxa(nano_bact_raref_2reads_phy_unsh,Phylum=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_phy_unsh %>% sample_sums() %>% sum())) # 100% of Nanopore reads assigned for unshared Phyla
1 - ((subset_taxa(nano_bact_raref_2reads_phy_sh,Phylum=="__")) %>% sample_sums() %>% sum() /
       (nano_bact_raref_2reads_phy_sh %>% sample_sums() %>% sum())) # 88.3% of Nanopore reads assigned for shared Phyla


illu_bact_raref_2reads_fam_unsh %>% sample_sums() %>% sum()
illu_bact_raref_2reads_ord_unsh %>% sample_sums() %>% sum()
illu_bact_raref_2reads_class_unsh %>% sample_sums() %>% sum()
illu_bact_raref_2reads_phy_unsh %>% sample_sums() %>% sum()





###################################@
### ARCHAEA
illu_arc_2reads <- microbiome::core(illu_arc, detection=1, prevalence=0) # remove singletons after rarefaction
nano_arc_raref <- phyloseq_coverage_raref(nano_arc, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across bact samples (0.9411947)
nano_arc_raref_2reads <- microbiome::core(nano_arc_raref, detection=1, prevalence=0) # remove singletons after rarefaction

illu_arc_2reads_phy <- tax_glom(illu_arc_2reads, taxrank = "Phylum")
phyla_illu <- illu_arc_2reads_phy@tax_table[,2]
illu_arc_2reads_class <- tax_glom(illu_arc_2reads, taxrank = "Class")
class_illu <- paste(illu_arc_2reads_class@tax_table[,2],illu_arc_2reads_class@tax_table[,3])
illu_arc_2reads_ord <- tax_glom(illu_arc_2reads, taxrank = "Order")
orders_illu <- paste(illu_arc_2reads_ord@tax_table[,2],illu_arc_2reads_ord@tax_table[,3],illu_arc_2reads_ord@tax_table[,4])
illu_arc_2reads_fam <- tax_glom(illu_arc_2reads, taxrank = "Family")
fam_illu <- paste(illu_arc_2reads_fam@tax_table[,2],illu_arc_2reads_fam@tax_table[,3],illu_arc_2reads_fam@tax_table[,4],
                  illu_arc_2reads_fam@tax_table[,5])
illu_arc_2reads_gen <- tax_glom(illu_arc_2reads, taxrank = "Genus")
genus_illu <- paste(illu_arc_2reads_gen@tax_table[,2],illu_arc_2reads_gen@tax_table[,3],illu_arc_2reads_gen@tax_table[,4],
                    illu_arc_2reads_gen@tax_table[,5],illu_arc_2reads_gen@tax_table[,6])
species_illu <- paste(illu_arc_2reads@tax_table[,2],illu_arc_2reads@tax_table[,3],illu_arc_2reads@tax_table[,4],
                      illu_arc_2reads@tax_table[,5],illu_arc_2reads@tax_table[,6],illu_arc_2reads@tax_table[,7])

nano_arc_raref_2reads_phy <- tax_glom(nano_arc_raref_2reads, taxrank = "Phylum")
phyla_nano <- nano_arc_raref_2reads_phy@tax_table[,2]
nano_arc_raref_2reads_class <- tax_glom(nano_arc_raref_2reads, taxrank = "Class")
class_nano <- paste(nano_arc_raref_2reads_class@tax_table[,2],nano_arc_raref_2reads_class@tax_table[,3])
nano_arc_raref_2reads_ord <- tax_glom(nano_arc_raref_2reads, taxrank = "Order")
orders_nano <- paste(nano_arc_raref_2reads_ord@tax_table[,2],nano_arc_raref_2reads_ord@tax_table[,3],nano_arc_raref_2reads_ord@tax_table[,4])
nano_arc_raref_2reads_fam <- tax_glom(nano_arc_raref_2reads, taxrank = "Family")
fam_nano <- paste(nano_arc_raref_2reads_fam@tax_table[,2],nano_arc_raref_2reads_fam@tax_table[,3],nano_arc_raref_2reads_fam@tax_table[,4],
                  nano_arc_raref_2reads_fam@tax_table[,5])
nano_arc_raref_2reads_gen <- tax_glom(nano_arc_raref_2reads, taxrank = "Genus")
genus_nano <- paste(nano_arc_raref_2reads_gen@tax_table[,2],nano_arc_raref_2reads_gen@tax_table[,3],nano_arc_raref_2reads_gen@tax_table[,4],
                    nano_arc_raref_2reads_gen@tax_table[,5],nano_arc_raref_2reads_gen@tax_table[,6])
species_nano <- paste(nano_arc_raref_2reads@tax_table[,2],nano_arc_raref_2reads@tax_table[,3],nano_arc_raref_2reads@tax_table[,4],
                      nano_arc_raref_2reads@tax_table[,5],nano_arc_raref_2reads@tax_table[,6],nano_arc_raref_2reads@tax_table[,7])


nano_arc_raref_2reads_phy@tax_table[,2] %>% length() # nb phyla detected by nano
nano_arc_raref_2reads_class@tax_table[,3] %>% length() # nb classes detected by nano
nano_arc_raref_2reads_ord@tax_table[,4] %>% length() # nb orders detected by nano
nano_arc_raref_2reads_fam@tax_table[,5] %>% length() # nb families detected by nano
nano_arc_raref_2reads_gen@tax_table[,6] %>% length() # nb genus detected by nano
nano_arc_raref_2reads@tax_table[,7] %>% length() # nb species detected by nano

table(str_count(nano_arc_raref_2reads@tax_table[,7], "__"))[1]/(nano_arc_raref_2reads@tax_table[,7] %>% length())
table(str_count(nano_arc_raref_2reads_gen@tax_table[,6], "__"))[1]/(nano_arc_raref_2reads_gen@tax_table[,6] %>% length())
table(str_count(nano_arc_raref_2reads_fam@tax_table[,5], "__"))[1]/(nano_arc_raref_2reads_fam@tax_table[,5] %>% length())
table(str_count(nano_arc_raref_2reads_ord@tax_table[,4], "__"))[1]/(nano_arc_raref_2reads_ord@tax_table[,4] %>% length())
table(str_count(nano_arc_raref_2reads_class@tax_table[,3], "__"))[1]/(nano_arc_raref_2reads_class@tax_table[,3] %>% length())
table(str_count(nano_arc_raref_2reads_phy@tax_table[,2], "__"))[1]/(nano_arc_raref_2reads_phy@tax_table[,2] %>% length())

illu_arc_2reads_phy@tax_table[,2] %>% length() # nb phyla detected by illu
illu_arc_2reads_class@tax_table[,3] %>% length() # nb classes detected by illu
illu_arc_2reads_ord@tax_table[,4] %>% length() # nb orders detected by illu
illu_arc_2reads_fam@tax_table[,5] %>% length() # nb families detected by illu
illu_arc_2reads_gen@tax_table[,6] %>% length() # nb genus detected by illu
illu_arc_2reads@tax_table[,7] %>% length() # nb species detected by illu

table(str_count(illu_arc_2reads@tax_table[,7], "__"))[1]/(illu_arc_2reads@tax_table[,7] %>% length())
table(str_count(illu_arc_2reads_gen@tax_table[,6], "__"))[1]/(illu_arc_2reads_gen@tax_table[,6] %>% length())
table(str_count(illu_arc_2reads_fam@tax_table[,5], "__"))[1]/(illu_arc_2reads_fam@tax_table[,5] %>% length())
table(str_count(illu_arc_2reads_ord@tax_table[,4], "__"))[1]/(illu_arc_2reads_ord@tax_table[,4] %>% length())
table(str_count(illu_arc_2reads_class@tax_table[,3], "__"))[1]/(illu_arc_2reads_class@tax_table[,3] %>% length())
table(str_count(illu_arc_2reads_phy@tax_table[,2], "__"))[1]/(illu_arc_2reads_phy@tax_table[,2] %>% length())

((illu_arc_2reads_phy@tax_table[,2] %>% length()) - setdiff(phyla_illu,phyla_nano) %>% length()) # nb phyla detected by illu only
((nano_arc_raref_2reads_phy@tax_table[,2] %>% length()) - setdiff(phyla_nano,phyla_illu) %>% length()) # nb phyla detected by nano only
((illu_arc_2reads_class@tax_table[,3] %>% length()) - setdiff(class_illu,class_nano) %>% length()) # nb classes detected by illu only
((nano_arc_raref_2reads_class@tax_table[,3] %>% length()) - setdiff(class_nano,class_illu) %>% length()) # nb classes detected by nano only
((illu_arc_2reads_ord@tax_table[,4] %>% length()) - setdiff(orders_illu,orders_nano) %>% length()) # nb orders detected by illu only
((nano_arc_raref_2reads_ord@tax_table[,4] %>% length()) - setdiff(orders_nano,orders_illu) %>% length()) # nb orders detected by nano only
((illu_arc_2reads_fam@tax_table[,5] %>% length()) - setdiff(fam_illu,fam_nano) %>% length()) # nb families detected by illu only
((nano_arc_raref_2reads_fam@tax_table[,5] %>% length()) - setdiff(fam_nano,fam_illu) %>% length()) # nb families detected by nano only
((illu_arc_2reads_gen@tax_table[,6] %>% length()) - setdiff(genus_illu,genus_nano) %>% length()) # nb genus detected by illu only
((nano_arc_raref_2reads_gen@tax_table[,6] %>% length())  - setdiff(genus_nano,genus_illu) %>% length())# nb genus detected by nano only
((illu_arc_2reads@tax_table[,7] %>% length()) - setdiff(species_illu,species_nano) %>% length()) # nb species detected by illu only
((nano_arc_raref_2reads@tax_table[,7] %>% length()) - setdiff(species_nano,species_illu) %>% length()) # nb species detected by nano only

illu_arc_2reads %>% sample_sums %>% max
nano_arc_raref_2reads %>% sample_sums %>% max

13/14
14/15
24/31
6/11
7/23
10/34
13/50
14/83
24/171

43/54
99/140
194/316
285/483
393/785
522/1495



#######################################################@
## Phylogenetic stats by sequencer


stats <- read.csv("/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/tables/stats_seq_taxa.csv") %>% as.data.frame()
colnames(stats) <- c("Rank", "Nb_taxa","seq","detected_by_the_seq_only","percent_of_taxa_dectected_by_seq","percent_of_taxa_shared_with_the_other_seq")
positions <- c("Phylum","Class","Order","Family","Genus","Species")

##pdf("/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/figures/phylogenetic_stats_by_sequencer.pdf")
ggplot(stats, aes(x=Rank, y=Nb_taxa)) + geom_bar(aes(fill=factor(seq)), stat="identity", position="dodge") + scale_x_discrete(limits = positions)
ggplot(stats, aes(x=Rank, y=detected_by_the_seq_only)) + geom_bar(aes(fill=factor(seq)), stat="identity", position="dodge") + 
  scale_x_discrete(limits = positions) + ylim(0,1500)
ggplot(stats, aes(x=Rank, y=percent_of_taxa_shared_with_the_other_seq, color=factor(seq), group=factor(seq))) + geom_line() + 
  geom_point() + scale_x_discrete(limits = positions) + ylim(0,100)
##dev.off()


stats_reads <- read.csv("/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/tables/stats_seq_reads.csv", check.names = TRUE, sep="\t") %>% as.data.frame()
colnames(stats_reads) <- c("Rank", "Nb_reads","seq","detected_by_the_seq_only","percent_of_reads_dectected_by_seq","percent_of_reads_shared_with_the_other_seq")
positions <- c("Phylum","Class","Order","Family","Genus","Species")

##pdf("/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/figures/phylogenetic_stats_reads_by_sequencer.pdf")
ggplot(stats_reads, aes(x=Rank, y=detected_by_the_seq_only)) + geom_bar(aes(fill=factor(seq)), stat="identity", position="dodge") + 
  scale_x_discrete(limits = positions)
ggplot(stats_reads, aes(x=Rank, y=percent_of_reads_shared_with_the_other_seq, color=factor(seq), group=factor(seq))) + geom_line() + 
  geom_point() + scale_x_discrete(limits = positions) + ylim(0,100)
##dev.off()

#######################################################@
## Venn diagrams
## install.packages("VennDiagram")
library("VennDiagram")


venn_phy <- draw.pairwise.venn(area1 = 45, area2 = 54, cross.area = 43, category=c("Illumina","Nanopore"), scaled=TRUE, 
                               cex=2, col=c('yellow', 'turquoise'), fill=c('orange', 'turquoise'))
venn_class <- draw.pairwise.venn(area1 = 106, area2 = 140, cross.area = 99, category=c("Illumina","Nanopore"), scaled=TRUE, 
                                 cex=2, col=c('yellow', 'turquoise'), fill=c('orange', 'turquoise'))
venn_ord <- draw.pairwise.venn(area1 = 209, area2 = 316, cross.area = 194, category=c("Illumina","Nanopore"), scaled=TRUE, 
                               cex=2, col=c('yellow', 'turquoise'), fill=c('orange', 'turquoise'))
venn_fam <- draw.pairwise.venn(area1 = 309, area2 = 483, cross.area = 285, category=c("Illumina","Nanopore"), scaled=TRUE, 
                               cex=2, col=c('yellow', 'turquoise'), fill=c('orange', 'turquoise'))
venn_gen <- draw.pairwise.venn(area1 = 448, area2 = 785, cross.area = 393, category=c("Illumina","Nanopore"), scaled=TRUE, 
                               cex=2, col=c('yellow', 'turquoise'), fill=c('orange', 'turquoise'))
venn_esp <- draw.pairwise.venn(area1 = 749, area2 = 1495, cross.area = 522, category=c("Illumina","Nanopore"), scaled=TRUE, 
                               cex=2, col=c('yellow', 'turquoise'), fill=c('orange', 'turquoise'))

##pdf("/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/figures/venn_diagrams_phylogenet.pdf", width = 12, height = 2)
##grid.newpage()
grid.arrange(venn_phy, venn_class, venn_ord, venn_fam, venn_gen, venn_esp, nrow=1)
##dev.off()


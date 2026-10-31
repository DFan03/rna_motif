# QC report
- raw dir: data/raw/orthoflavi_v1
- records read: 238
- records kept: 207
- removed/flagged by reason: {'FLAG_moltype_DNA_cDNA': 10, 'non_ACGU_frac>0.01': 13, 'length_outlier_MADz>3.5': 16, 'FLAG_no_CDS_annotation': 12, 'UTR_incomplete': 1, 'exact_duplicate_of': 1}
- params: max_amb_frac=0.01, len_mad_z=3.5, min_scale=50.0
- versions: python 3.12.7, biopython 1.84, pandas 2.2.2

## Kept per species
species_query
Tick-borne encephalitis virus       20
Murray Valley encephalitis virus    19
West Nile virus                     19
Yellow fever virus                  19
Dengue virus 2                      18
Dengue virus 4                      18
Zika virus                          18
Dengue virus 3                      17
Japanese encephalitis virus         15
Saint Louis encephalitis virus      15
Usutu virus                         15
Dengue virus 1                      14

## Length (nt) per species
                                  count      min      50%      max
species_query                                                     
Dengue virus 1                     14.0  10511.0  10680.0  10745.0
Dengue virus 2                     18.0  10477.0  10678.5  10724.0
Dengue virus 3                     17.0  10640.0  10694.0  10709.0
Dengue virus 4                     18.0  10593.0  10644.5  10653.0
Japanese encephalitis virus        15.0  10871.0  10965.0  10977.0
Murray Valley encephalitis virus   19.0  10865.0  10948.0  11014.0
Saint Louis encephalitis virus     15.0  10762.0  10929.0  10940.0
Tick-borne encephalitis virus      20.0  10531.0  10899.5  11104.0
Usutu virus                        15.0  10728.0  10954.0  11066.0
West Nile virus                    19.0  10174.0  10910.0  11029.0
Yellow fever virus                 19.0  10395.0  10862.0  11007.0
Zika virus                         18.0  10636.0  10764.5  10808.0

## UTR lengths (nt)
         utr5_len    utr3_len
count  195.000000  195.000000
mean    91.830769  480.251282
std     21.686308  117.907282
min     24.000000   74.000000
25%     79.000000  412.500000
50%     95.000000  451.000000
75%    101.000000  570.000000
max    142.000000  727.000000

## Missing metadata (kept, flagged)
missing_host                42
missing_collection_date     26
missing_country            207
missing_taxid                0

## Hosts
host
Homo sapiens               91
NA                         42
Culex annulirostris        12
Turdus merula               6
Culex tarsalis              5
mosquito                    3
Ixodes persulcatus          3
Culex tritaeniorhynchus     2

table MSGbed
"Multisector Glioma Illumina EPIC methylation bed file format"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"Name of gene"
uint    score;		"Beta value"
uint  	reserved;	"Green indicates unmethylated, red indicates methylated"
)


## ---- createGR ----
library(GenomicRanges)
gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),
           strand=c("+","-","-")
)
gr
# subset like a data frame
gr[1:2,]


## ---- createGRwMetadata ----
gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),
           names=c("id1","id3","id2"),
           scores=c(100,90,50)
)
# or add it later (replaces the existing meta data)
mcols(gr)=DataFrame(name2=c("pax6","meis1","zic4"),
                    score2=c(1,2,3))

gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),
           names=c("id1","id3","id2"),
           scores=c(100,90,50)
)

# or appends to existing meta data
mcols(gr)=cbind(mcols(gr),
                          DataFrame(name2=c("pax6","meis1","zic4")) )
gr
# elementMetadata() and values() do the same things
elementMetadata(gr)
values(gr)

## ---- convertDataframe2gr ----
# read CpGi data set
cpgi.df = read.table("data/cpgi.hg19.chr21.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
# remove chr names with "_"
cpgi.df =cpgi.df [grep("_",cpgi.df[,1],invert=TRUE),]

cpgi.gr=GRanges(seqnames=cpgi.df[,1],
                ranges=IRanges(start=cpgi.df[,2],
                              end=cpgi.df[,3]))
# read refseq file
ref.df = read.table("data/refseq.hg19.chr21.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
ref.gr=GRanges(seqnames=ref.df[,1],
               ranges=IRanges(start=ref.df[,2],
                              end=ref.df[,3]),
               strand=ref.df[,6],name=ref.df[,4])
# get TSS
tss.gr=ref.gr
# end of the + strand genes must be equalized to start pos
end(tss.gr[strand(tss.gr)=="+",])  =start(tss.gr[strand(tss.gr)=="+",])
# startof the - strand genes must be equalized to end pos
start(tss.gr[strand(tss.gr)=="-",])=end(tss.gr[strand(tss.gr)=="-",])
# remove duplicated TSSes ie alternative transcripts
# this keeps the first instance and removes duplicates
tss.gr=tss.gr[!duplicated(tss.gr),]

## ---- findTSSwithCpGi ----
subsetByOverlaps(tss.gr,cpgi.gr)


## ---- countOverlaps ----
counts=countOverlaps(tss.gr,cpgi.gr)  


## ---- findOverlaps ----
findOverlaps(tss.gr,cpgi.gr)  

## ---- findNearest ----
# find nearest CpGi to each TSS
n.ind=nearest(tss.gr,cpgi.gr)
# get distance to nearest
dists=distanceToNearest(tss.gr,cpgi.gr,select="arbitrary")
dists

# histogram of the distances to nearest TSS
dist2plot=mcols(dists)[,1]
hist(log10(dist2plot),xlab="log10(dist to nearest TSS)",
     main="Distances")


## ---- findCanonical ----
  


## ---- cannonicalCoverage ----


## ---- getCoverageBam ----
# regions of interest
# promoters on chr21
promoter.gr=tss.gr
start(promoter.gr)=start(promoter.gr)-1000
end(promoter.gr)  =end(promoter.gr)+1000
promoter.gr=promoter.gr[seqnames(promoter.gr)=="chr21"]

library(Rsamtools)
bamfile="data/wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam"

# get reads for regions of interest from the bam file
param <- ScanBamParam(which=promoter.gr)
alns <- readGAlignmentsFromBam(bamfile, param=param)

covs=coverage(alns) # get coverage vectors

myViews=Views(covs,as(promoter.gr,"RangesList")) # get subsets of coverage
cov.mat=t(viewApply(myViews$chr21,as.vector)) # get coverage matrix

# plot mean coverage
plot(-1000:1000,colMeans(cov.mat),type="l",
     ylab="mean coverage",xlab="window around TSS",
     main="ChIP-Seq SP1")

#image(cov.mat)

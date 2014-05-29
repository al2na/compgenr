Title
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **Help** toolbar button for more details on using R Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:



```r
library(GenomicRanges)
gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),
           strand=c("+","-","-")
)
gr
```

```
## GRanges with 3 ranges and 0 metadata columns:
##       seqnames     ranges strand
##          <Rle>  <IRanges>  <Rle>
##   [1]     chr1 [ 50, 100]      +
##   [2]     chr2 [150, 200]      -
##   [3]     chr2 [200, 300]      -
##   ---
##   seqlengths:
##    chr1 chr2
##      NA   NA
```

```r
# subset like a data frame
gr[1:2,]
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames     ranges strand
##          <Rle>  <IRanges>  <Rle>
##   [1]     chr1 [ 50, 100]      +
##   [2]     chr2 [150, 200]      -
##   ---
##   seqlengths:
##    chr1 chr2
##      NA   NA
```


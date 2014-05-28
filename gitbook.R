
#devtools::install_github("al2na/Rgitbook")
require(knitr)
require(Rgitbook)


#buildRmd(dir = getwd(), clean = FALSE )

# my (Altuna Akalin) on buildGitbook
# since I want to edit in gitbook app not Rstudio, I keep
# the summary files have .Rmd locations not .md
# this function temporarily replaces SUMMARY.md
# and then runs gitbook
buildGitbook2<-function(source.dir = getwd(), out.dir = paste0(getwd(), "/_book"),
                        buildRmd = TRUE, format, gitbook.params,...){
  
  sums=readLines(paste0(source.dir,"/SUMMARY.md"))
  file.rename(paste0(source.dir,"/SUMMARY.md"), paste0(getwd(),"/.SUMMARY.md"))
  message('Writing new SUMMARY.md...')
  f <- file(paste0(source.dir,"/SUMMARY.md"))
  writeLines(gsub(".Rmd",".md",sums), f)
  close(f)
  
  buildGitbook(source.dir = source.dir, out.dir = out.dir,
               buildRmd = buildRmd, format, gitbook.params, ...)
  file.rename(paste0(source.dir,"/.SUMMARY.md"), paste0(source.dir,"/SUMMARY.md"))
  
}


buildGitbook2(source.dir = getwd(), out.dir = paste0(getwd(), "/_book"),
              buildRmd = TRUE )

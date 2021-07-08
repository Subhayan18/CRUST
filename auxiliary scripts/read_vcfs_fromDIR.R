directory <- "<directory name here>"

temp <-   list.files(path = directory, pattern = "*.vcf", 
                     full.names = T)
myfiles = lapply(temp, read.vcfR)

names<-gsub(".*/(.+).vcf.*", "\\1", temp)
for (i in 1 : length(names)){ assign(names[i],myfiles[[i]]) }
rm('temp','myfiles','i','names')

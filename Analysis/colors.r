colors <- list()
colors[['Genotype']] <- c(RColorBrewer::brewer.pal(n=3, 'Set1')[1:2])
names(colors[['Genotype']]) <- c('E3', 'E4')

colors[['Sex']] <- c(RColorBrewer::brewer.pal(n=3, "Set2")[1:2])
names(colors[['Sex']]) <- c('M', 'F')

colors[['Group']] <- c(RColorBrewer::brewer.pal(n=4, 'Set3'))
names(colors[['Group']]) <- c('E3_F', 'E3_M', 'E4_F', 'E4_M')

colors[['lipid.class']] <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n=8, 'Dark2'))(13)
names(colors[['lipid.class']]) <- c("CL", "LCL", "PS", "PG", "SM", "PC", "LPC", "PE", "LPE", "CAR", "FFA", "CE", "free")
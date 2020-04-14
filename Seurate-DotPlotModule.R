
DotPlotModule=function(object, assay = NULL, geneList, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
          group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, scale.max = NA)
{
  assay <- ifelse(is.null(assay),DefaultAssay(object = object),assay)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.geneList <- lapply(geneList, function(x) FetchData(object = object, vars = x))
  names(data.geneList)=names(geneList)
  data.geneList$id <- if (is.null(x = group.by)) { Idents(object = object) 
    } else { object[[group.by, drop = TRUE]] }
  if (!is.factor(x = data.geneList$id)) {
    data.geneList$id <- factor(x = data.geneList$id)
  }
  id.levels <- levels(x = data.geneList$id)
  data.geneList$id <- as.vector(x = data.geneList$id)
  data.plot <- lapply(X = unique(x = data.geneList$id), FUN = function(ident) {
    selected=data.geneList$id == ident
    data.use = lapply(data.geneList[1:length(data.geneList)-1], function(x){
      x[selected,]
    })
    avg.exp <- sapply(X = data.use, FUN = function(x) {
        return(mean(x = rowMeans(expm1(x = x))))
      # return(mean(x = expm1(x = rowMeans(x))))
      })
    threshold=0
    pct.exp <- sapply(X = data.use, function(x){
      x=as.numeric(rowMeans(x))
      return(length(x = x[x > threshold]) / length(x = x))
    })
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.geneList$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min, 
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = names(geneList)))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
                                      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
                         2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
                     no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

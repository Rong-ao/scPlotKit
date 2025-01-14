#' @title Orthogonal DotPlot
#' @description Drawing DotPlot showing gene expression levels across two classification groups
#' @param object seurat object
#' @param features gene list to show
#' @param assay seurat object assay to draw, defaults to the active assay
#' @param cols customized color list
#' @param col.min Minimum scaled average expression threshold
#' @param col.max Maximum scaled average expression threshold
#' @param dot.min set minimum size of dots
#' @param dot.scale The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param idents Identity classes to include in plot (default is all)
#' @param cluster.idents Whether to order identities by hierarchical clusters
#' based on given features, default is FALSE
#' @param scale Determine whether the data is scaled, TRUE for default
#' @param legend Combine legends into a single legend
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param keep.scale How to handle the color and size scale across multiple plots. Same as Seurat::FeaturePlot()
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#' @param by.col If splitting by a factor, plot the splits per column with the features as rows
#' @param ncols For several features, set number of columns in merged plot
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' @param group.by.x First factor to group the cells by
#' @param group.by.y Second factor to group the cells by
#'
#' @return A ggplot object
#'
#' @export
#' @concept visualization
#'
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom cowplot theme_cowplot
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size scale_radius
#' theme element_blank labs scale_color_identity scale_color_distiller element_rect
#' dup_axis guides element_text margin scale_color_brewer scale_color_gradientn
#' scale_color_manual coord_fixed ggtitle scale_color_gradient guides
#' guide_legend guide_colorbar scale_x_continuous scale_y_continuous theme
#' facet_grid unit
#' @importFrom scattermore geom_scattermore
#' @importFrom stats dist hclust
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom Seurat DefaultAssay CellsByIdentities FetchData
#' @importFrom dplyr
#'
#' @export
#' @concept visualization
#'
#' @aliases SplitDotPlotGG
#' @seealso \code{RColorBrewer::brewer.pal.info}
#'
#' @examples
#' ifnb <- LoadData("ifnb")
#' OrthoDotPlot(ifnb, features = c("CD4", "CD86", "CD96"), group.by.x = "stim", group.by.y = "seurat_annotations", legend = T, keep.scale = 'feature')
#' OrthoDotPlot(ifnb, features = c("CD4", "CD86", "CD96"), group.by.x = "stim", group.by.y = "seurat_annotations", legend = T, keep.scale = 'all')
OrthoDotPlot <- function(
    object,
    features,
    assay = NULL,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    legend = NULL,
    scale.by = 'radius',
    keep.scale = 'feature',
    scale.min = NA,
    scale.max = NA,
    coord.fixed = FALSE,
    by.col = TRUE,
    ncols = NULL,
    combine = TRUE,
    group.by.x,
    group.by.y){
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  # split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  split.colors <- F
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), idents = idents))
  data.features <- FetchData(object = object, vars = features, cells = cells)
  group.by <- group.by.x
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  data.features$x_id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  data.features$x_id <- factor(x = data.features$x_id)
  data.features$y_id <- if (is.null(x = group.by.y)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by.y, drop = TRUE]][cells, drop = TRUE]
  }
  data.features$y_id <- factor(x = data.features$y_id)

  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = group.by.y)) {
    splits <- FetchData(object = object, vars = group.by.y)[cells, group.by.y]
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }

  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      m = 2
      if (!is.null(x = group.by.y)) {
        m = 3
      }
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - m), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          # Average expression in Dotplot was calculated as mean of
          # "exp(x)-1 transformed" expression of all cells (raw counts)
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      data.use$x_id <- unique(subset(data.features, id == x)$x_id)
      data.use$y_id <- unique(subset(data.features, id == x)$y_id)
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = log1p(data.use))
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  color.by <- 'avg.exp.scaled'
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }

  plots <- vector(mode = "list", length = length(x = features))
  for (j in 1:length(x = features)) {
    feature <- features[j]
    plot <- ggplot(data = subset(data.plot, features.plot == feature), mapping = aes_string(x = 'x_id', y = 'y_id')) +
      geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
      scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(size = guide_legend(title = 'Percent Expressed')) +
      labs(
        x = 'Group.X',
        y = 'Group.Y'
      ) +
      theme_cowplot()
    if (!is.null(x = feature.groups)) {
      plot <- plot + facet_grid(
        facets = ~feature.groups,
        scales = "free_x",
        space = "free_x",
        switch = "y"
      ) + theme(
        panel.spacing = unit(x = 1, units = "lines"),
        strip.background = element_blank()
      )
    }
    if (length(x = cols) == 1) {
      plot <- plot + scale_color_distiller(palette = cols)
    } else {
      plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
    plot <- plot + labs(title = feature)
    if (!(is.null(x = keep.scale)) && keep.scale == "feature") {
      max.feature.value <- max(subset(data.plot, features.plot == feature)$avg.exp.scaled)
      min.feature.value <- min(subset(data.plot, features.plot == feature)$avg.exp.scaled)
      plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
    }
    # Add coord_fixed
    if (coord.fixed) {
      plot <- plot + coord_fixed()
    }
    plot <- plot
    plots[[j]] <- plot
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  # Combine the plots
  if (is.null(x = ncols)) {
    ncols <- 2
    if (length(x = features) == 1) {
      ncols <- 1
    }
    if (length(x = features) > 6) {
      ncols <- 3
    }
    if (length(x = features) > 9) {
      ncols <- 4
    }
  }
  ncols <- ifelse(
    test = is.null(x = group.by.y),
    no = ncols,
    yes = length(x = features)
  )
  # legend <- group.by.y %iff% 'none'
  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (isTRUE(x = combine)) {
    plots <- wrap_plots(plots, ncols = ncols, nrow = group.by.y %iff% 1)
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all") {
      max.feature.value <- max(data.plot$avg.exp.scaled)
      min.feature.value <- min(data.plot$avg.exp.scaled)
      max.percent.value <- max(data.plot$pct.exp)
      min.percent.value <- min(data.plot$pct.exp)
      # plots <- plots & NoLegend()
      plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)) &
                                  scale.func(range = c(0, dot.scale), limits = c(min.percent.value, max.percent.value)))
      plots <- plots + plot_layout(guides='collect')
    }
  }
  return(plots)
}

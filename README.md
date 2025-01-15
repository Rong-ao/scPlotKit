# scPlotKit
This R package is the collection of some customized little tools for scRNA data visualization.
Functions here are mostly based on source code of Seurat v5, see https://github.com/satijalab/seurat
More functions will come soon...

Thanks for my dummy friend, Jiayi Song in Xiamen University, for his endless requirments of visualization in scRNA-seq data analysis pushing me develop this tiny package.

***2025.1.14: First submission of scPlotKit***

This kit now contains one visualization function: `OrthoDotPlot()`, which is constrcuted based on `Seurat::DotPlot()`. This function provides a direct way to draw a gene expression level DotPlot with two classifications (metadata columns) of cells in seurat object, showing each as an axis of DotPlot and corresponding gene expression pattern in orthogonal cell subgroup. It also supports to draw several genes in batch and generate combined plot with separated scale or not.

**Install**: 
`devtools::install_github("Rong-ao/scPlotKit")`
or
`remotes::install_github("Rong-ao/scPlotKit")`

**`OrthoDotPlot()` Usage example:**

```
library(scPlotKit)
libaray(Seurat)
library(SeuratData)
ifnb <- LoadData("ifnb")
p1 <- OrthoDotPlot(ifnb, features = c("CD4", "CD86", "CD96"),
                   group.by.x = "stim", group.by.y = "seurat_annotations",
                   legend = T, keep.scale = 'feature')
p2 <- OrthoDotPlot(ifnb, features = c("CD4", "CD86", "CD96"),
                   group.by.x = "stim", group.by.y = "seurat_annotations",
                   legend = T, keep.scale = 'all')
```

Plots below show p1 and p2 from code above:
![image](https://github.com/user-attachments/assets/dfd6163f-4c1d-40ce-85d4-77a64494c596) ![image](https://github.com/user-attachments/assets/34abc107-203b-4c4d-b173-a3fa0f9f9f2b)

You can check more arguments of `OrthoDotPlot()` in R for more nice plots!

If any issue, please contact with kourongao@westlake.edu.cn
